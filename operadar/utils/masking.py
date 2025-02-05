#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from operadar.utils.make_links import link_keys_with_available_hydrometeors


def mask_precipitations(contents:dict[np.ndarray],expMmin:float,hydrometeorMoments:dict[int]):
    "Mask grid points where the total precipitations are below expMmin threshold."
    print("Masking precipitations where total <",10**expMmin)
    hydrometeors_model = link_keys_with_available_hydrometeors(hydrometeorMoments,datatype='model')
    total_content=np.sum(contents[h] for h in hydrometeors_model)
    mask_precip=(total_content>10**expMmin)  
    return mask_precip   
