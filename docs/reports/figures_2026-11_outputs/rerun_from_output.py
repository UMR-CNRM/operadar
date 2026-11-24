"""Rerun operadar (current code) from the model fields stored in a previous operadar output, when the AROME .fa
file cannot be read (no FA library). Only the reading of the model file is replaced; tables, geometry, mixed
phase, dual-pol variables, sensitivity threshold and saving are the operadar code. Attenuation (Zh_att) is not
computed (pressure and humidity are not stored in the outputs).

    cd <repo> && python docs/reports/figures_2026-11_outputs/rerun_from_output.py OLD_OUTPUT.nc BAND OUT_DIR
"""
import sys
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, '.')
import operadar.forward_operator as fo
from operadar.load_config_file import load_configuration_file

src, band, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
ds = xr.open_dataset(src)
assert ds.attrs['mixed_phase_type'] == 'Fw_posg', 'contents before the mixed phase are only recoverable with Fw_posg'


def read_model_file(filePath, modelname, micro_scheme, real_case, domain, hydrometeorMoments, verbose=False):
    C = {str(h): ds.Contents.sel(hydrometeor=h).values.astype('f8') / 1000. for h in ds.hydrometeor.values}  # kg/m3
    N = {str(h): ds.Concentrations.sel(hydrometeor=h).values.astype('f8') for h in ds.hydrometeor.values}
    # Fw_posg: in the melting layer the graupel was moved into wg (gg = 0 there), elsewhere wg = 0
    C['gg'] = C['gg'] + C.pop('wg')
    N.pop('wg')
    M = {h: C[h] for h in ['cc', 'rr', 'ss', 'gg', 'ii'] if h in C}
    Nc = {h: N[h] for h in M}
    X, Y = ds.x.values.astype('i4'), ds.y.values.astype('i4')
    return [X, Y, ds.Alt.values.astype('f8'), ds.lon.values, ds.lat.values, M, Nc, ds['T'].values.astype('f8'),
            None, None]


fo.read_model_file = read_model_file
fo.format_temporal_variable = lambda filePath, model_type, real_case: pd.Timestamp(ds.time.values)

fo.operadar(filename='historic.arome.franmg-01km30+0018:00.fa', configuration=load_configuration_file('test_AROME'),
            radar_band=band, tables_path='./tables_generator/tables/David2026PhD/',
            out_dir_path=out_dir, compute_attenuation=False, subDomain=None)
