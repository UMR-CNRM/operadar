#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Statistics behind the CFADs (Contoured Frequency by Altitude Diagrams, Yuter and Houze 1995).

compute_cfad.py stores, for each variable and each height, the histogram of the values with fine bins
(VAR_INFO[...]['hist'] in operadar_output.py). Everything is derived from these counts by plot_cfad.py:
- quantiles (median, Q25-Q75...) for the CFADs of David et al. (2025, AMT), with a precision of one fine bin;
- frequency CFADs (2D histograms) with coarser bins (sums of fine bins).
Counts of several files (cases, periods) can simply be added.
"""
import numpy as np
import xarray as xr

from .operadar_output import var_info


def fine_edges(name: str) -> np.ndarray:
    """Fine bins of the histograms: (min, max, width) for regular bins, ('log', first, max, n) for the contents
    (0, then n bins regularly spaced in log between `first` and `max`)."""
    h = var_info(name)['hist']
    if h[0] == 'log':
        return np.concatenate([[0.], np.logspace(np.log10(h[1]), np.log10(h[2]), h[3])])
    lo, hi, w = h
    n = int(round((hi - lo) / w))
    return np.round(lo + w * np.arange(n + 1), 10)  # exact edges (0 is an edge)


class CfadAccumulator:
    """Histograms (height x bin) of several variables, filled file after file."""

    def __init__(self, heights, names):
        self.heights = np.asarray(heights, 'f4')
        self.edges = {n: fine_edges(n) for n in names}
        self.counts = {n: np.zeros((len(self.heights), len(e) - 1), np.int64) for n, e in self.edges.items()}
        # exact zeros are counted apart (e.g. no rain above the melting layer, Zdr = 0 of spheres), so that a
        # median of zeros is 0 and not somewhere inside the first bin
        self.zeros = {n: np.zeros(len(self.heights), np.int64) for n in self.edges}

    def add(self, name: str, values: np.ndarray):
        """values: (n_heights, ...) with NaN where the point is not selected. Values outside the histogram
        range are counted in the first or last bin (so that the quantiles stay right)."""
        e = self.edges[name]
        nb = len(e) - 1
        values = values.reshape(len(self.heights), -1)
        for i in range(len(self.heights)):
            v = values[i]
            v = v[np.isfinite(v)]
            is_zero = v == 0
            if is_zero.any() and e[0] <= 0 <= e[-1]:
                self.zeros[name][i] += int(is_zero.sum())
                v = v[~is_zero]
            if v.size:
                idx = np.clip(np.searchsorted(e, v, side='right') - 1, 0, nb - 1)
                self.counts[name][i] += np.bincount(idx, minlength=nb)

    def to_dataset(self, attrs: dict) -> xr.Dataset:
        ds = xr.Dataset(coords={'height': ('height', self.heights, {'units': 'm'})}, attrs=attrs)
        for n in self.counts:
            ds[f'edges_{n}'] = (f'edge_{n}', self.edges[n].astype('f8'))
            ds[f'count_{n}'] = (('height', f'bin_{n}'), self.counts[n],
                                {'long_name': f'number of values of {n} in each bin (exact zeros excluded)'})
            ds[f'zeros_{n}'] = ('height', self.zeros[n], {'long_name': f'number of values of {n} equal to 0'})
        return ds


def names_in(stats: xr.Dataset) -> list[str]:
    return [v[len('count_'):] for v in stats.data_vars if v.startswith('count_')]


def _zero_bin(edges: np.ndarray) -> int:
    return int(np.clip(np.searchsorted(edges, 0., side='right') - 1, 0, len(edges) - 2))


def get(stats: xr.Dataset, name: str):
    """(counts (height, bin), edges) of a variable, exact zeros put back in the bin of 0."""
    counts, edges = stats[f'count_{name}'].values.copy(), stats[f'edges_{name}'].values
    if f'zeros_{name}' in stats:
        counts[:, _zero_bin(edges)] += stats[f'zeros_{name}'].values
    return counts, edges


def get_for_quantiles(stats: xr.Dataset, name: str):
    """(counts, edges) with the exact zeros as a bin of zero width [0, 0], for quantiles()."""
    counts, edges = stats[f'count_{name}'].values, stats[f'edges_{name}'].values
    if f'zeros_{name}' not in stats:
        return counts, edges
    k = _zero_bin(edges)
    e = np.concatenate([edges[:k], [0., 0.], edges[k + 1:]])
    c = np.concatenate([counts[:, :k], stats[f'zeros_{name}'].values[:, None], counts[:, k:]], axis=1)
    return c, e


def merge(stats_list: list[xr.Dataset]) -> xr.Dataset:
    """Sum the counts of several statistics files (same heights and bins)."""
    out = stats_list[0].copy(deep=True)
    for s in stats_list[1:]:
        if not np.array_equal(s.height.values, out.height.values):
            raise ValueError('cannot merge statistics computed on different heights')
        for n in names_in(s):
            for v in (f'count_{n}', f'zeros_{n}'):
                if v in out and v in s:
                    out[v] = out[v] + s[v]
                elif v in s:
                    out[v] = s[v]
            if f'edges_{n}' not in out:
                out[f'edges_{n}'] = s[f'edges_{n}']
            elif not np.allclose(out[f'edges_{n}'].values, s[f'edges_{n}'].values):
                raise ValueError(f'cannot merge {n}: different bins')
        for key in ('n_files', 'n_times', 'n_mask_columns'):
            if key in s.attrs and key in out.attrs:
                out.attrs[key] = int(out.attrs[key]) + int(s.attrs[key])
    return out


def quantiles(counts: np.ndarray, edges: np.ndarray, percents, min_count: int = 1) -> np.ndarray:
    """Percentiles (n_heights, n_percents) from histograms, linear inside the bin; NaN when fewer than
    min_count values at that height."""
    percents = np.atleast_1d(percents)
    out = np.full((counts.shape[0], len(percents)), np.nan)
    for i, c in enumerate(counts):
        total = c.sum()
        if total < max(min_count, 1):
            continue
        cum = np.cumsum(c)
        for j, p in enumerate(percents):
            target = p / 100. * total
            k = int(np.searchsorted(cum, target, side='left'))
            k = min(k, len(c) - 1)
            before = cum[k - 1] if k > 0 else 0
            frac = (target - before) / c[k] if c[k] > 0 else 0.
            out[i, j] = edges[k] + frac * (edges[k + 1] - edges[k])
    return out


def rebin(counts: np.ndarray, edges: np.ndarray, width: float, lim=None, log=False):
    """Sum the fine bins into coarser bins of `width` (or `width` decades if log), between lim=(min, max)
    (default: the range of the fine bins). Each fine bin goes to the coarse bin containing its centre;
    fine bins outside lim are dropped (they stay in the totals used for the normalisation per level)."""
    centres = 0.5 * (edges[:-1] + edges[1:])
    if log:
        lo, hi = lim if lim is not None else (edges[1], edges[-1])
        coarse = 10 ** np.arange(np.log10(lo), np.log10(hi) + width / 2, width)
    else:
        lo, hi = lim if lim is not None else (edges[0], edges[-1])
        lo = np.floor(lo / width + 1e-9) * width
        coarse = np.round(np.arange(lo, hi + width / 2, width), 10)
    j = np.searchsorted(coarse, centres, side='right') - 1
    ok = (j >= 0) & (j < len(coarse) - 1)
    assign = np.zeros((len(centres), len(coarse) - 1))
    assign[np.nonzero(ok)[0], j[ok]] = 1.
    return counts @ assign, coarse
