#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 13:17:47 2018
@author: augrosc & lebastardt & montangonm & sinhorin & davidcl
"""

# External modules
import argparse
import time as tm
from pathlib import Path

# 0perad modules
import operadar.operadar_conf as cf
from operadar.read.model import read_model_file
from operadar.utils.masking import mask_precipitations
from operadar.radar.geometry import compute_radar_geometry
from operadar.read.tmatrix_tables import read_Tmatrix_Clotilde
from operadar.microphysics.mixed_phase import compute_mixed_phase
from operadar.radar.dualpol_variables import compute_dualpol_variables
from operadar.utils.make_links import link_keys_with_available_hydrometeors
from operadar.save.save_dpolvar import create_tree_structure_outFiles, save_netcdf
from operadar.utils.formats_data import format_temporal_variable,define_output_path



def operadar(filename:str,
             read_tmatrix:bool=True,
             out_dir_path:str=cf.outPath,
             radar_band:str=cf.radar_band,
             Tmatrix_params:dict={},
             subDomain:list[float]|None=cf.subDomain,
             get_more_details=False,
             ) -> tuple[bool,bool,dict]:
    """Radar forward operator main function. 

    Args:
        filename (str): Only the name of the file
        read_tmatrix (bool, optional): Option to save computing time within a loop. Need to be set to True at least for the first iteration. Defaults to True.
        out_dir_path (str, optional): Path to store the output files. Defaults to cf.outPath.
        radar_band (str, optional): Defaults to cf.radar_band.
        Tmatrix_params (dict, optional): Dictionnary containing the Tmatrix tables parameters. Defaults to {}.
        subDomain (list[float] | None, optional): Defaults to cf.subDomain.
        get_more_details (bool): Defaults to False

    Returns:
        read_tmatrix (bool): to update the boolean during a loop
        Tmatrix_params (dict) : to keep in memory the Tmatrix parameters and values throughout multiple iterations over the same radar band
        
    """
    
    begin_program_timer = tm.time()
    
    # Create or check tree structure of the output directory path
    create_tree_structure_outFiles(output_dir=Path(out_dir_path))
    
    # Format temporal variable and output file name
    input_file_path = Path(cf.input_filePath+filename)
    temporal_variable = format_temporal_variable(filePath=input_file_path)
    outFilePath = define_output_path(out_dir_path,radar_band,temporal_variable) 
    
    if not Path(outFilePath).with_suffix('.nc').exists():
        
        # Read Tmatrix tables (files from Clotilde)
        if read_tmatrix :
            Tmatrix_hydromet_list = link_keys_with_available_hydrometeors(hydrometeorMoments=cf.hydrometeors_moments,
                                                                          datatype='tmatrix',
                                                                          )
            Tmatrix_params = read_Tmatrix_Clotilde(band=radar_band,
                                                   hydrometeors=Tmatrix_hydromet_list,
                                                   verbose=get_more_details,
                                                   )
        # Read model variables
        [X, Y, Z, lon, lat, M, Nc, Tc] = read_model_file(filePath=input_file_path,
                                                         domain=subDomain,
                                                         verbose=get_more_details,
                                                         )
        # Compute radar geometry
        mask_dist_max, elevations = compute_radar_geometry(X=X,
                                                           Y=Y,
                                                           Z=Z,
                                                           Tc=Tc,
                                                           elev_max=Tmatrix_params['ELEVmax']["rr"])
        # Mask precipitations
        mask_precip = mask_precipitations(contents=M,
                                          expMmin=Tmatrix_params['expMmin']["rr"],
                                          hydrometeorMoments=cf.hydrometeors_moments)
        # Combine masks
        partial_mask = (mask_precip & mask_dist_max)
        
        # Compute mixed phase parametrization
        [M, Nc,Fw] = compute_mixed_phase(contents=M,
                                         concentrations=Nc,
                                         expMmin=Tmatrix_params['expMmin']["rr"]) 
        # Compute dual-pol radar variables
        dpolDict = compute_dualpol_variables(temperature=Tc,
                                             mask_precip_dist=partial_mask,
                                             elev=elevations, Fw=Fw,
                                             contents=M,
                                             concentrations=Nc,
                                             tmatrix_param=Tmatrix_params,
                                             output_file_path=outFilePath,
                                             X=X, Y=Y, Z=Z,
                                             lat=lat, lon=lon,
                                             date_time=temporal_variable)
        # Saving file
        save_netcdf(X=X, Y=Y, Z=Z, lat=lat, lon=lon,
                    datetime=temporal_variable, dpolDict=dpolDict,
                    contentsDict=M, concentrationsDict=Nc,
                    temperature=Tc, outfile=Path(outFilePath),
                    )
        # For multiple iterations over different time but with the same settings, save time by not reading
        # again Tmatrix tables and lat/lon fields (if available)
        read_tmatrix = False
        elapsed_time = tm.time() - begin_program_timer
        print("Elapsed time :",int(elapsed_time//60),"minutes",int(elapsed_time%60),"seconds")
        print("-----------------------------------------------------------------------------")
        return read_tmatrix, Tmatrix_params
    
    else :
        print("File exists at :",outFilePath+'.nc')
        print("-----------------------------------------------------------------------------")
        read_tmatrix = True
        return read_tmatrix, dict()




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="OPERADAR")
    parser.add_argument("filename", type=str)
    parser.add_argument("--verbose", action="store_true",default=False)
    operad_args = parser.parse_args()
    operadar(filename=operad_args.filename, get_more_details=operad_args.verbose)
