#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import datetime as dt
from operadar.forward_operator import operadar
from operadar.load_config_file import load_configuration_file


# stratiforme case with 0Z run and convective case with 12Z run
cases = [('conv','2022-06-20 03:00','2022-06-20 11:00')] # ('strat','2024-11-19 06:00','2024-11-19 22:00'),
microphysics_schemes = ['ICJW','LIMC'] # 'ICE3','ICJW','LIMASG','LIMC'

for case_type,begin,end in cases :
    for micro in microphysics_schemes :
        # Read tables and loop over the datetime (again) when changing microphysics and case stuy     
        read_tables = True 
        dict_tables = {}
        ech = pd.to_datetime(begin, format="%Y-%m-%d %H:%M")
        end = pd.to_datetime(end, format="%Y-%m-%d %H:%M")
        # Reload the corresponding configuration in the python environment
        config = load_configuration_file(f'conf_{micro}_{case_type}.py')
        # Execute operadar for each arome file
        while ech <= end :
            time = ech.strftime('%H:%M')
            print('\n','------------------------------------------',case_type,micro,time,'------------------------------------------')
            fname = f'historic.arome.franmg-01km30+00{time}.fa'
            read_tables, dict_tables = operadar(filename=fname,
                                                configuration=config,
                                                read_tables=read_tables,
                                                tables_content = dict_tables,
                                                get_more_details=False,
                                                )
            ech += dt.timedelta(minutes=5)

