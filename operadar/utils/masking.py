#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import operadar.operad_conf as cf
from operadar.utils.make_links import link_keys_with_available_hydrometeors


def mask_precipitations(contents:dict[np.ndarray],
                        expMmin:float,
                        hydrometeorMoments:dict[int],
                        )-> np.ndarray:
    """Mask grid points where the total precipitations are below 10^expMmin."""
    print("Masking precipitations where total <",10**expMmin,"kg/m3.")
    hydrometeors_model = link_keys_with_available_hydrometeors(hydrometeorMoments=cf.hydrometeors_moments,datatype='model')
    total_content=np.sum(contents[h] for h in hydrometeors_model)
    mask_precip= total_content>10**expMmin  
    return mask_precip   



def mask_bright_band(contents:dict[np.ndarray],
                     expMmin:float,
                     )-> np.ndarray:
    """Mask grid points where rain, graupel and/or hail contents are below 10^expMmin."""
    hydrometeors_model = link_keys_with_available_hydrometeors(hydrometeorMoments=cf.hydrometeors_moments,datatype='model')
    if 'hh' in hydrometeors_model :
        mask_BB = ((contents["rr"] > 10**expMmin) & ((contents["gg"]> 10**expMmin) | (contents["hh"]> 10**expMmin)))
    else :							
        mask_BB=((contents["rr"] > 10**expMmin) & (contents["gg"]> 10**expMmin))
    return mask_BB



def mask_hydrometeor(content:np.ndarray,
                     expMmin:float,
                     )-> np.ndarray:
    """Mask grid points where hydrometeor contents are below 10^expMmin."""
    mask_hydro = content > 10**expMmin
    return mask_hydro