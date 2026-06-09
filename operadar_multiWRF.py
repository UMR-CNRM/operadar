#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from operadar.forward_operator import operadar
from operadar.load_config_file import load_configuration_file

import glob
import os
import xarray as xr
from pathlib import Path
import shutil 

#Check if there exist the WRF/tmp directory. If not, create it. If it exists, delete its contents.

dir_tmp = Path("modelFiles/WRF/tmp/") 
dir_tmp.mkdir(parents=True, exist_ok=True)
for filename in os.listdir(dir_tmp):
    file_path = os.path.join(dir_tmp, filename)
    try:
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)
    except Exception as e:
        print('Failed to delete %s. Reason: %s' % (file_path, e))


# stratiforme case with 0Z run and convective case with 12Z run
cases = ['241025'] # ('strat','2024-11-19 06:00','2024-11-19 22:00'),
microphysics_schemes = ['WRFICE3'] # 'ICE3','ICJW','LIMASG','LIMC'

files = sorted(glob.glob("modelFiles/WRF/split_files/wrfout_t*.nc"))
files = [os.path.basename(f) for f in files]
for case_type in cases :
    for micro in microphysics_schemes :
        config = load_configuration_file(f'conf_{micro}_{case_type}.py') 
        for file_name in files:

            read_tables = True
            dict_tables = {}

            print('processing...', file_name)
            read_tables, dict_tables = operadar(filename=file_name,
                                                configuration=config,
                                                read_tables=read_tables,
                                                tables_content = dict_tables,
                                                get_more_details=False,
                                                )
            
print('Operadar Finished. Joining files....')

all_files = set(glob.glob("modelFiles/WRF/split_files/dpolvar*.nc")) #The dpolvar files that I have now

new_files = sorted(list(all_files))
datasets = [xr.open_dataset(f) for f in new_files]
ds = xr.concat(datasets, dim='time')
ds.to_netcdf('modelFiles/WRF/split_files/dpolvar_WRF_ICE3_Wband_joined.nc')