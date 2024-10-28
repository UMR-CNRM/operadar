#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import datetime as dt
import operad_conf as cf


# ========== Read csv configuration file ==========
def read_csv_file(csv_file_path:str,
                  which_date:str,
                  time_columns_name:list=["start_time", "end_time"],
                  csv_datetime_format:str="%Y-%m-%d %H:%M",
                  csv_delimiter=";"
                  ):
    df = pd.read_csv(csv_file_path, delimiter=csv_delimiter)
    df[time_columns_name] = df[time_columns_name].apply(lambda x: pd.to_datetime(x, format=csv_datetime_format))
    if which_date != "all" :
        date  = dt.datetime.strptime(which_date, "%Y%m%d").date()
        subdf = df.loc[df[time_columns_name[0]].dt.date == date]
        return subdf
    else :
        return df


# ========== Extract domain given in the csv file ==========
def extract_csv_domain(csv_row,domain_columns_name):
    if type(domain_columns_name)==str:
        domain = [float(x) for x in csv_row[domain_columns_name].strip().split(',')]
        lon_min = domain[0] ; lon_max = domain[1]
        lat_min = domain[2] ; lat_max = domain[3]
    elif type(domain_columns_name)==list:       # rajouter une condition pour vérifier que les 4 sous éléments sont des str
        lon_min = csv_row[domain_columns_name[0]] ; lon_max = csv_row[domain_columns_name[1]]
        lat_min = csv_row[domain_columns_name[2]] ; lat_max = csv_row[domain_columns_name[3]]
    else :
        print('Unsupported domain format')
        #print(csv_instructions)
    return lat_min,lat_max,lon_min,lon_max


# ========== Extract global info (common to MNH and AROME) ==========
def extract_csv_info(csv_row,time_columns_name:list,model_run_column_name:str ,microphysics_scheme):
    deb = csv_row[time_columns_name[0]]
    fin = csv_row[time_columns_name[1]]
    radar_band = str(csv_row.radar_band)
    
    # --------- Cloe's special --------- # (no impact on the creation of the file)
    try :
        radar_ids = "-".join([str(x) for x in csv_row.radar_id_list.strip().split(',')])
    except :
        radar_ids = ""
    # ---------------------------------- #
    
    if model_run_column_name != None :
        run  = int(csv_row[model_run_column_name])
        output_path = f"{cf.outPath}/{deb.strftime('%Y%m%d')}/{str(run).zfill(2)}Z_{microphysics_scheme}_k{cf.MixedPhase}/{radar_ids}"
    else :
        output_path = f"{cf.outPath}/{deb.strftime('%Y%m%d')}/{microphysics_scheme}_k{cf.MixedPhase}/{radar_ids}"
    
    return deb,fin,run,radar_band,output_path
