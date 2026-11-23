<!-- To add to the wiki page "Simulation configuration" (section on the radar options) -->

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
