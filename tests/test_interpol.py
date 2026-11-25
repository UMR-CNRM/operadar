"""The hypercube interpolation of operadar must match an independent multilinear interpolation of the lookup
table (scipy RegularGridInterpolator on the axes Tc, ELEV, log10 N or Fw, log10 M), and give the table values
exactly on the grid nodes.

The former reference implementation (CALC_KTMAT + INTERPOL, kept in this file until November 2026) read the
row of the next content (index + 1): it is replaced by this independent check.
"""
from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import RegularGridInterpolator

from operadar.read.interpolate_tables import hypercube_interpolation
from operadar.read.lookup_tables import read_and_extract_tables_content

TABLES = Path(__file__).resolve().parents[1] / 'tables_generator' / 'tables' / 'David2026PhD'
BAND = 'C'
DPOL = ['Zh', 'Zdr', 'Kdp']
# (hydrometeor, scheme, moment, P3 name used by operadar)
CASES = [('ii', 'ICE3', 1, 'Nc'), ('wg', 'ICE3', 1, 'Fw'), ('rr', 'LIMA', 2, 'Nc')]


def load(h, scheme, moment):
    name = f'TmatCoefInt_{scheme}_{moment}M_{BAND}{h}'
    if not any((TABLES / f'{name}{ext}').exists() for ext in ('', '.nc')):
        pytest.skip(f'{name} not available in {TABLES}')
    tables, _, columns = read_and_extract_tables_content(band=BAND, hydrometeors=[h], moments={h: moment},
                                                         scheme=scheme, dpol2add=DPOL, path_table=f'{TABLES}/',
                                                         verbose=False, test_interpolation=True)
    tables['Nc'] = tables['N']
    # scattering coefficients only (the grid columns Tc, ELEV, M, Fw/N are not interpolated)
    columns = sorted(c for c in columns if c not in ('Tc', 'ELEV', 'M', 'P3', 'Fw', 'N')
                     and h in tables.get(c, {}))
    assert columns, 'no scattering coefficient read from the table'
    return tables, columns


def axes(t, h, p3name):
    def ax(lo, hi, step):
        return lo + step * np.arange(int(round((hi - lo) / step)) + 1)
    pk = 'Fw' if p3name == 'Fw' else 'expCC'
    return (ax(t['Tcmin'][h], t['Tcmax'][h], t['Tcstep'][h]), ax(t['ELEVmin'][h], t['ELEVmax'][h], t['ELEVstep'][h]),
            ax(t[f'{pk}min'][h], t[f'{pk}max'][h], t[f'{pk}step'][h]),
            ax(t['expMmin'][h], t['expMmax'][h], t['expMstep'][h]))


def random_points(grid, n, rng):
    return [rng.uniform(g[0], g[-1], n) if len(g) > 1 else np.full(n, g[0]) for g in grid]


@pytest.mark.parametrize('h,scheme,moment,p3name', CASES)
def test_hypercube_matches_independent_interpolation(h, scheme, moment, p3name):
    t, columns = load(h, scheme, moment)
    grid = axes(t, h, p3name)
    shape = tuple(len(g) for g in grid)
    rng = np.random.default_rng(0)
    T, E, P, lM = random_points(grid, 500, rng)
    P3 = P if p3name == 'Fw' else 10 ** P
    new = hypercube_interpolation(t, h, p3name, E.copy(), T.copy(), P3.copy(), 10 ** lM, columns)
    pts = np.column_stack([T, E, P, lM])
    keep = [k for k, g in enumerate(grid) if len(g) > 1]  # Fw has a single value for dry species
    for col in columns:
        values = np.asarray(t[col][h], dtype=float)
        assert values.size == np.prod(shape), f'{col}: table size {values.size} != grid {shape}'
        cube = values.reshape(shape)
        ref = RegularGridInterpolator([grid[k] for k in keep], cube.reshape([shape[k] for k in keep]))(pts[:, keep])
        np.testing.assert_allclose(new[col], ref, rtol=1e-6, atol=1e-30, err_msg=f'{h} {col}')


@pytest.mark.parametrize('h,scheme,moment,p3name', CASES)
def test_hypercube_returns_table_values_on_nodes(h, scheme, moment, p3name):
    t, columns = load(h, scheme, moment)
    grid = axes(t, h, p3name)
    shape = tuple(len(g) for g in grid)
    rng = np.random.default_rng(1)
    idx = [rng.integers(0, n, 200) for n in shape]
    T, E, P, lM = (g[i] for g, i in zip(grid, idx))
    P3 = P if p3name == 'Fw' else 10 ** P
    new = hypercube_interpolation(t, h, p3name, E.copy(), T.copy(), P3.copy(), 10 ** lM, columns)
    flat = np.ravel_multi_index(idx, shape)
    for col in columns:
        np.testing.assert_allclose(new[col], np.asarray(t[col][h], dtype=float)[flat], rtol=1e-6, atol=1e-30,
                                   err_msg=f'{h} {col}')
