#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import datetime as dt
from operadar.forward_operator import operadar
from operadar.load_config_file import load_configuration_file


# stratiforme case with 0Z run and convective case with 12Z run
cases = [('241025', '2024-10-25 12:00', '2024-10-26 23:00')] # ('strat','2024-11-19 06:00','2024-11-19 22:00'),
microphysics_schemes = ['ICE3'] # 'ICE3','ICJW','LIMASG','LIMC'

for case_type,begin,end in cases :
    for micro in microphysics_schemes :
        # Read tables and loop over the datetime (again) when changing microphysics and case study     

        ech = pd.to_datetime(begin, format="%Y-%m-%d %H:%M")
        end = pd.to_datetime(end, format="%Y-%m-%d %H:%M")
        # Reload the corresponding configuration in the python environment
        config = load_configuration_file('conf_AROICE3_241025.py')
        begin = pd.to_datetime(begin, format="%Y-%m-%d %H:%M")

        # Execute operadar for each arome file
        while ech <= end :
            read_tables = True 
            dict_tables = {}

            dT = int((ech-begin).total_seconds()/3600)
            time = ech.strftime('%H:%M')

            print('\n','------------------------------------------',case_type,micro,time,'------------------------------------------')
            fname = f'historic.arome.franmg-01km30+{dT:04d}:00.fa'
            read_tables, dict_tables = operadar(filename=str(fname),
                                                configuration=config,
                                                read_tables=read_tables,
                                                tables_content = dict_tables,
                                                get_more_details=False,
                                                )
            ech += dt.timedelta(minutes=60)

