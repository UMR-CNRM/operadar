#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import datetime as dt
from operadar.forward_operator import operadar
from operadar.load_config_file import load_configuration_file


# CASES FRANCE
# cases = [('FRANCE','2023-05-09','historic.arome.franmg-01km30+0000:00.fa'),
#          ('FRANCE','2023-06-08','historic.arome.franmg-01km30+0001:00.fa'),
#          ('FRANCE','2023-06-19','historic.arome.franmg-01km30+0004:00.fa'),
#          ('FRANCE','2023-06-22','historic.arome.franmg-01km30+0001:00.fa'),
#          ('FRANCE','2023-07-09','historic.arome.franmg-01km30+0004:00.fa'),]

# CASES ANTILLES
# cases = [('ANTILLES','2023-06-25','historic.arome.caraibes-01km30+0002:00.fa'),
#          ('ANTILLES','2023-06-29','historic.arome.caraibes-01km30+0002:00.fa'),
#          ('ANTILLES','2023-08-21','historic.arome.caraibes-01km30+0002:00.fa'),          
#          ('ANTILLES','2023-09-09','historic.arome.caraibes-01km30+0004:00.fa'),
#          ('ANTILLES','2023-10-20','historic.arome.caraibes-01km30+0010:00.fa'),
#          ('ANTILLES','2024-11-08','historic.arome.caraibes-01km30+0003:00.fa'),
cases = [('ANTILLES','2023-10-20','historic.arome.caraibes-01km30+0010:00.fa'),
         ('ANTILLES','2024-11-08','historic.arome.caraibes-01km30+0003:00.fa'),
         ('CALEDONIA','2023-03-23','historic.arome.caledonie-01km30+0005:00.fa'),
         ('CALEDONIA','2023-05-30','historic.arome.caledonie-01km30+0002:00.fa'),
         ('CALEDONIA','2025-01-06','historic.arome.caledonie-01km30+0001:00.fa'),         
         ('POLYNESIA','2023-05-11','historic.arome.polynesie-01km30+0004:00.fa'),
         ('POLYNESIA','2023-09-18','historic.arome.polynesie-01km30+0005:00.fa'), 
         ('POLYNESIA','2024-03-15','historic.arome.polynesie-01km30+0004:00.fa'), 
         ('POLYNESIA','2024-04-21','historic.arome.polynesie-01km30+0004:00.fa'), 
         ('POLYNESIA','2025-01-03','historic.arome.polynesie-01km30+0005:00.fa'),         
         ('REUNION','2023-04-11','historic.arome.indien-01km30+0000:00.fa'), 
         ('REUNION','2023-12-19','historic.arome.indien-01km30+0001:00.fa'), 
         ('REUNION','2023-12-23','historic.arome.indien-01km30+0005:00.fa'), 
         ('REUNION','2024-01-13','historic.arome.indien-01km30+0003:00.fa'), 
         ('REUNION','2024-01-16','historic.arome.indien-01km30+0001:00.fa'), 
         ('REUNION','2024-01-17','historic.arome.indien-01km30+0004:00.fa'), 
         ('REUNION','2024-05-29','historic.arome.indien-01km30+0002:00.fa'), 
         ('REUNION','2024-06-12','historic.arome.indien-01km30+0003:00.fa'), 
         ]
         

TmatConfig='David2026PhD' # 'default'
configFile=f'test_AROME_Antia_{TmatConfig}.py'


for domain,day,fname in cases :   
    #read_tables = True # REMOVE THIS LINE IF THE TABLES RED ARE TO BE THE SAME FOR EACH CASE
    modelFilesDir=f'/cnrm/obs/data1/raspaudd/GNSS-PRO/ANTIA/{domain}/{day}/AROME/'
    outputDir=f'/home/augros/Programmes/operadar/operadarFiles/{TmatConfig}/{domain}/'
    dict_tables = {}
    config = load_configuration_file(configFile) # PROVIDE HERE A CONFIG FILE COMMON TO ALL CASES (DO NOT CARE ABOUT THE subDomain ARGUMENT IN THE CONFIG FILE, IT WILL BE OVERWRITTEN


    print('\n','------------------------------------------',domain,day,fname,'------------------------------------------')
    read_tables, dict_tables = operadar(filename=fname,
                                        configuration=config,
                                        in_dir_path=modelFilesDir,
                                        out_dir_path=outputDir,
                                        tables_content = dict_tables,
                                        get_more_details=False,
                                        )
