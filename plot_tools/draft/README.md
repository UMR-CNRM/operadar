# Former scripts (kept for reference)

These scripts were written for specific cases (paths of C. David and C. Augros, CORSE case) and some need
libraries that are not in the repository (`lib/cfad_lib.py`, `lib/level2alt.py`). Their features are now in the
scripts of `plot_tools/`, which read the operadar outputs directly (see `plot_tools/README.md`):

| Draft script | Replaced by |
|---|---|
| `compute_CFADs.py`, `exec_computeCFAD.sh`, `stats_for_CFADs_cloe.py` (quantiles per altitude inside a mask, pickle files) | `compute_cfad.py` (statistics files in NetCDF, any mask, many files) |
| `plot_CFADs.py`, `tracer_CFADs.ipynb` (median and Q25-Q75 of the radar variables and contents) | `plot_cfad.py --figs summary radar contents fw` |
| `cfad_frequency.py`, `cfad_frequency2.py` (frequency CFADs, differences with the observations) | `plot_cfad.py --figs freq --diff [--log]` |
| `plot_2Dcumulmaps.py` (maximum over several times with the 40 dBZ contours) | `plot_max_over_time.py` |
| `plot_D_M_varpol_all_tempe.py` (curves of one type at several temperatures) | `plot_tables_2D.py` (M-Tc maps) |
| `plot_D_M_varpol_Rayleigh_test.py` | `plot_D_M_varpol_Rayleigh.py` |

Differences with the drafts: interpolation of the model levels on heights inside the scripts (no
`*interp.nc` files), cell cores and envelopes computed from the column maximum of Zh as in David et al. (2025,
AMT), quantiles from histograms (constant memory, statistics of several cases can be added), Fw of the wet
graupel rebuilt from the contents.
