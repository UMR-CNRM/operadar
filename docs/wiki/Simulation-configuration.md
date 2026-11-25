# Configuration file
The user must set up, beforehand, a configuration file based on the template provided `./configFiles/template.py`. Please, do not modify the template directly. Instead, make a copy of the file in `./configFiles/` and name it differently. Examples ready to use are `test_AROME_Cband.py`, `test_AROME_Wband_vertical.py` and `test_MesoNH_Cband.py` (see [Execution of operadar](https://github.com/UMR-CNRM/operadar/wiki/Execution-of-operadar)).

The options of the configuration file are:

| Parameter      | Description |
| -----------    | ----------- |
| `input_filePath` | Path to the folder containing the model file(s), e.g. `"./modelFiles/AROME/"` |
| `output_filePath` | Folder where the output files are stored, e.g. `"./operadarFiles/AROME/"` |
| `path_tables` | Lookup tables directory, e.g. `"./tables_generator/tables/David2026PhD/"` (see [Lookup table configurations](https://github.com/UMR-CNRM/operadar/wiki/Lookup-table-configurations)) |
| `model` | Model name: `'Arome'` or `'MesoNH'` |
| `real_case` | `True` or `False` (MesoNH idealized cases) |
| `microphysics_scheme` | Microphysics scheme: `'ICE3'`, `'ICJW'`, `'LIMA'` or `'LIMC'` + an optional name extension (e.g. `'LIMA_noHail'` or `'ICE3_CIBU_moins'`); only the four first characters select the tables |
| `hydrometeors_moments` | Number of moments of each hydrometeor of the scheme, e.g. `{'cc':1,'rr':1,'ss':1,'gg':1,'ii':1,'wg':1}` (add `'hh':1,'wh':1` for hail and wet hail, L to X bands) |
| `cloud_water_over` | `'land'` or `'sea'`: the cloud water tables `cl` (land) or `cs` (sea) are used |
| `subDomain` | Written as one of the following : <ul><li>`[lon_min,lon_max,lat_min,lat_max]` if `real_case=True`</li> <li>`[i_min,i_max,j_min,j_max]` if `real_case=False`</li> <li>`None` (all the points in the file). Large MesoNH domains may not fit in memory: use a subdomain</li> </ul>|
| `mixed_phase_parametrization` | How the mixed phase is simulated: `'T_pos'` (melting species only at positive temperatures), `'Fw_pos'` (rain and graupel moved into wet graupel in the melting layer) or `'Fw_posg'` (only the graupel moved into wet graupel) |
| `save_netcdf_single_hydrometeor` | If `True`, also compute the dual-pol variables of each hydrometeor type and save them (1 file per type) |
| `dpol2add` | List of at least one variable among `['Zh','Zdr','Kdp','Rhohv','Ah','Av']` |
| `scattering_method` | `'Tmatrix'` (only method currently used) |
| `compute_attenuation` | If `True`, operadar computes the attenuated reflectivity `'Zh_att'` (dBZ) of a vertically pointing radar. Needs `'Ah'` in `dpol2add` |
| `radar_altitude` | Altitude of the vertically pointing radar (m): 0 for a radar at the ground, or the aircraft altitude for an airborne radar |
| `sensitivity_threshold` | If `True` (default), the variables are removed where Zh is below the minimum detectable reflectivity (see [Radar sensitivity](#radar-sensitivity)) |
| `Zmin_1km` | Minimum detectable reflectivity at 1 km (dBZ); `None` uses a default per band |
| `radar_band` | Radar band name (see [Available radar bands](https://github.com/UMR-CNRM/operadar/wiki#available-radar-bands)) |
| `distmax_rad` | Maximum range (m) of the simulated radar, used when `radarloc` is set |
| `radarloc` | Radar location (MesoNH only; not used for AROME):<ul><li>`'center'` places the radar at the center of the domain</li><li>`[lat_radar,lon_radar]` places the radar at the given coordinates</li><li>`None` does not simulate the radar beam geometry</li></ul> |
| `cnst_angle` | Used when `radarloc=None` to simulate either a horizontal (`cnst_angle=0`) or a vertically pointing (`cnst_angle=90`) radar at each grid point |

### Attenuation geometry

The attenuation (`compute_attenuation = True`) is integrated **along the vertical** from `radar_altitude`, whatever
`cnst_angle`: it is the attenuation seen by a vertically pointing radar (at the ground looking up, or airborne).
Use `cnst_angle = 90` so that Zh, Zdr and Kdp are consistent with that geometry (with `cnst_angle = 0`, Zdr and Kdp
are those of a horizontal beam). Example: `configFiles/test_AROME_Wband_vertical.py`.

### Radar sensitivity

```python
sensitivity_threshold = True
Zmin_1km = None
```

A real radar does not measure anything where the reflectivity is below its minimum detectable reflectivity, and
the polarimetric variables computed from an extremely small signal (rhohv especially) are meaningless. With
`sensitivity_threshold = True`, operadar sets **Zh, Zdr, Kdp and Rhohv** (and `Zh_att`) to NaN where Zh is below

    Zmin(r) = Zmin_1km + 20 log10(r / 1 km)     [dBZ]

where r is the distance to the radar: it is known when the radar geometry is simulated (`radarloc`), and for a
vertically pointing radar (`cnst_angle = 90`: r = |altitude - radar_altitude|, zenith or nadir). Otherwise
(horizontal beam at every grid point) the threshold is `Zmin_1km` everywhere. When the attenuation is computed,
the attenuated reflectivity `Zh_att` is compared with Zmin. Ah and Av are kept.

`Zmin_1km = None` uses a default per band (dBZ at 1 km), an order of magnitude for typical radars:

| Band | L, S, C, X | Ku, K, Ka | W |
|---|---|---|---|
| Zmin at 1 km (dBZ) | -30 | -40 | -50 |

Set `Zmin_1km` to the value of the simulated radar (e.g. from its minimum detectable reflectivity at a given
range r0: Zmin_1km = Zmin(r0) - 20 log10(r0 / 1 km)). `sensitivity_threshold = False` keeps all the values.
