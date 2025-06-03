#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time as tm
import numpy as np
from operadar.utils.masking import mask_bright_band
from operadar.utils.make_links import link_keys_with_available_hydrometeors



def compute_mixed_phase(contents:dict[np.ndarray],
                        concentrations:dict[np.ndarray],
                        hydrometeorMoments:dict[int],
                        expMmin:float,
                        parametrization:str,
                        ) -> tuple[dict[np.ndarray],dict[np.ndarray],np.ndarray]:
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
    
    mask_BB = mask_bright_band(contents, hydrometeorMoments, expMmin)
    
    hydrometeors = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments,
                                                         datatype='tables',
                                                         )
    wet_species = [key for key in hydrometeors if (key[0:1]=="w") and (key[1:2]*2 in hydrometeors)]
    
    Fw = compute_liquid_water_fraction(contents,mask_BB,wet_species)
    
    for spec in wet_species :
        contents[spec] = np.copy(contents[spec[1:2]*2])
        concentrations[spec] = np.copy(concentrations[spec[1:2]*2])
    
    contents = apply_mixed_phase_parametrization(contents, Fw, mask_BB, wet_species,parametrization) 
    
    print("\t--> Done in",round(tm.time()- deb_timer,2),"seconds")       
    return contents, concentrations, Fw



def compute_liquid_water_fraction(contents:dict[np.ndarray],
                                  mask_BB:np.ndarray,
                                  wet_species:list
                                  ) -> np.ndarray:
    """Computing the liquid water fraction depending on the wet species in the config file."""
    print("\tComputing the liquid water fraction for wet species :",wet_species)
    Fw = np.zeros(np.shape(contents["rr"]))
    Fw[mask_BB] = (contents["rr"]/ (contents["rr"]+np.sum(contents[key[1:2]*2] for key in wet_species) ))[mask_BB]
    return Fw



def apply_mixed_phase_parametrization(contents:dict[np.ndarray],
                                      Fw:np.ndarray,
                                      BB:np.ndarray,
                                      wet_species:list,
                                      parametrization:str,
                                      ) -> dict[np.ndarray]:
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
     
    if parametrization == "Fw_pos" :
        contents["wg"][Fw == 0] = 0
        contents["wg"][BB] = contents["gg"][BB]+contents["rr"][BB] # If contents["rr"] > 10**expMmin) & (contents["gg"]> 10**expMmin)
                                   # the rainwater is added to the wet graupel content         
        contents["rr"][BB] = 0        # and removed from the rain content  
        contents["gg"][BB] = 0
       
        if ('wh' in wet_species) and ('hh' in wet_species) :
            contents["wh"] = contents["hh"] + ( (contents["rr"]*contents["hh"])/(contents["hh"]+contents["gg"]) )	# d'après Wolfensberger, 2018

    if parametrization == "Fw_posg" :
        contents["wg"][Fw == 0] = 0
        contents["wg"][BB] = contents["gg"][BB]
        contents["gg"][BB] = 0
       
        if ('wh' in wet_species) and ('hh' in wet_species) :
            contents["wh"][Fw == 0] = 0
            contents["wh"][BB] = contents["hh"][BB]
            contents["hh"][BB] = 0 # addition Clotilde (in the mixed phase, there is no dry hail)
    
    return contents