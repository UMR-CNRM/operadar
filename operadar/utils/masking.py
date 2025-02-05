#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import operadar.operad_conf as cf
from operadar.utils.make_links import link_keys_with_available_hydrometeors


def mask_precipitations(contents:dict[np.ndarray], expMmin:float, hydrometeorMoments:dict[int])-> np.ndarray:
    """Mask grid points where the total precipitations are below expMmin threshold."""
    print("Masking precipitations where total <",10**expMmin,"kg/m3.")
    hydrometeors_model = link_keys_with_available_hydrometeors(hydrometeorMoments=cf.moments,datatype='model')
    total_content=np.sum(contents[h] for h in hydrometeors_model)
    mask_precip=(total_content>10**expMmin)  
    return mask_precip   



def mask_bright_band(contents:dict[np.ndarray], expMmin:float)-> np.ndarray:
    """Mask grid points where rain, graupel and/or hail contents are below expMmin threshold."""
    hydrometeors_model = link_keys_with_available_hydrometeors(hydrometeorMoments=cf.moments,datatype='model')
    if 'hh' in hydrometeors_model :
        mask_BB = ((contents["rr"] > 10**expMmin) & ((contents["gg"]> 10**expMmin) | (contents["hh"]> 10**expMmin)))
    else :							
        mask_BB=((contents["rr"] > 10**expMmin) & (contents["gg"]> 10**expMmin))
    return mask_BB