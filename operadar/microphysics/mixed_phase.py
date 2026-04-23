#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time as tm
import numpy as np
from operadar.utils.masking import mask_bright_band
from operadar.utils.make_links import link_keys_with_available_hydrometeors



def compute_mixed_phase(contents:dict[str,np.ndarray],
                        concentrations:dict[str,np.ndarray],
                        hydrometeorMoments:dict[str,int],
                        expMmin:float,
                        parametrization:str,
                        ) -> tuple[dict[str,np.ndarray],dict[str,np.ndarray],np.ndarray]:
    """Compute the mixed phase following the parametrization set in the configuration file.

    Args:
        contents (dict[np.ndarray]): dict of 3D arrays (one per hydrometeor).
        concentrations (dict[np.ndarray]): dict of 3D arrays (one per hydrometeor).
        hydrometeorMoments (dict of form {str : int})): dict containing the number of moments for each hydrometeor of the microphysics scheme.
        expMmin (float): value stored in the table dict.
        parametrization (str): Mixed phase parametrization (either T_pos, Fw_pos or Fw_posg).

    Returns:
        contents (dict[np.ndarray]): dict of 3D arrays.
        concentrations (dict[np.ndarray]): dict of 3D arrays.
        Fw (np.ndarray): liquid water fraction (3D) .
    """
    
    print('Estimating mixed phase where rain water coexists with iced species (even at negative temperatures).') ; deb_timer = tm.time()
    
    mask_BB = mask_bright_band(contents=contents,
                               hydrometeors_moments=hydrometeorMoments,
                               expMmin=expMmin,
                               )
    
    hydrometeors = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments,
                                                         datatype='tables',
                                                         )
   # More robust way to find the dry counterpart (e.g., 'wg' -> 'gg')
    wet_species = [k for k in hydrometeors if k.startswith('w') and (k[1] * 2)
                   in hydrometeors]
    
    Fw = compute_liquid_water_fraction(contents, mask_BB, wet_species)
    
    for wspec in wet_species:
        dry_spec = wspec[1] * 2  # Converts 'wg' to 'gg'
        contents[wspec] = np.copy(contents[dry_spec])
        concentrations[wspec] = np.copy(concentrations[dry_spec])
    
    
    contents = apply_mixed_phase_parametrization(contents, Fw, mask_BB, wet_species,parametrization) 
    
    print("\t--> Done in",round(tm.time()- deb_timer,2),"seconds")       
    return contents, concentrations, Fw



def compute_liquid_water_fraction(contents:dict[str,np.ndarray],
                                  mask_BB:np.ndarray,
                                  wet_species:list
                                  ) -> np.ndarray:
    """Computing the liquid water fraction depending on the wet species in the config file."""
    print("\tComputing the liquid water fraction for wet species :",wet_species)

    Fw = np.zeros_like(contents["rr"])
    
    # Sum the dry species element-wise
    sum_dry = np.zeros_like(contents["rr"])
    for wspec in wet_species:
        sum_dry += contents[wspec[1] * 2]
    
    # Prevent division by zero by adding a tiny epsilon or using a mask
    denominator = contents["rr"] + sum_dry
    denominator[denominator == 0] = 1e-12 
    
    Fw[mask_BB] = contents["rr"][mask_BB] / denominator[mask_BB]
    return Fw



def apply_mixed_phase_parametrization(contents:dict[str,np.ndarray],
                                      Fw:np.ndarray,
                                      BB:np.ndarray,
                                      wet_species:list,
                                      parametrization:str,
                                      ) -> dict[str,np.ndarray]:
    """Available parametrizations :
    * T_pos : the species content is transferred to the melting species only at positive temperatures.
    * Fw_pos : the rain and graupel content are emptied and transferred into the wet graupel content within the melting layer.
    * Fw_posg : only the graupel content is emptied and transferred to the wet graupel content within the melting layer.
    """
    
    #if parametrization == "T_pos" :  
    #    contents["wg"][Tc < 0] = 0
    #    contents["gg"][Tc >= 0] = 0
    #    contents["wh"][Tc < 0] = 0
    #    contents["hh"][Tc >= 0] = 0
     
    if parametrization == "Fw_pos":
        contents["wg"][Fw == 0] = 0
        # 1. Calculate the new wet species FIRST, while 'rr' and 'gg' still have their values
        if ('wh' in wet_species) and ('hh' in wet_species):
            denom = contents["hh"][BB] + contents["gg"][BB]
            denom[denom == 0] = 1e-12 # Safety
            contents["wh"][BB] = contents["hh"][BB] + ((contents["rr"][BB] * contents["hh"][BB]) / denom[BB])
        
        # 2. Now perform the transfers/zeroing
        contents["wg"][BB] = contents["gg"][BB] + contents["rr"][BB]
        contents["rr"][BB] = 0        
        contents["gg"][BB] = 0
        
    if parametrization == "Fw_posg" :
        contents["wg"][Fw == 0] = 0
        contents["wg"][BB] = contents["gg"][BB]
        contents["gg"][BB] = 0
       
        if ('wh' in wet_species) and ('hh' in wet_species) :
            contents["wh"][Fw == 0] = 0
            contents["wh"][BB] = contents["hh"][BB]
            contents["hh"][BB] = 0 # addition Clotilde (in the mixed phase, there is no dry hail)
    
    return contents