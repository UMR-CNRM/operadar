#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import epygram
import time as tm
import numpy as np
from pathlib import Path

epygram.init_env()



def append_in_input_file(complete_input_path:Path,
                         dpolVar:dict[str,np.ndarray],
                         var2add:list,
                         ):
    
    loaded_epygram_file = epygram.formats.resource(filename=complete_input_path,
                                                   openmode = 'a',
                                                   fmt = 'FA',
                                                   )
    levels = loaded_epygram_file.geometry.vcoordinate.levels
    deb = tm.time()
    for lvl in levels:
        tmp_field = loaded_epygram_file.readfield(f'S0{str(lvl).zfill(2)}RAIN')
        for var in var2add:
            field_copy=tmp_field.clone() 
            field_copy.fid['FA']=f'S0{str(lvl).zfill(2)}_{var}' # rename
            data2replace =  field_copy.data*0 # set everything to 0
            data2replace += dpolVar[var][lvl-1]
            field_copy.setdata(data2replace)
            loaded_epygram_file.writefield(field_copy)
    print('Done in',np.round(tm.time()-deb,2),'seconds')
    loaded_epygram_file.close()