#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 13:17:47 2018
@author: augrosc & lebastardt & montangonm & sinhorin & davidcl
"""

# External modules
import argparse
import time as tm
from sys import exit
from pathlib import Path

# 0perad modules
import operadar.operadar_conf as cf
from operadar.read.model import read_model_file
from operadar.utils.masking import mask_precipitations
from operadar.radar.geometry import compute_radar_geometry
from operadar.save.append_dpolvar import append_in_input_file
from operadar.read.tmatrix_tables import read_Tmatrix_Clotilde
from operadar.microphysics.mixed_phase import compute_mixed_phase
from operadar.radar.dualpol_variables import compute_dualpol_variables
from operadar.utils.make_links import link_keys_with_available_hydrometeors
from operadar.save.save_dpolvar import create_tree_structure_outFiles, save_netcdf
from operadar.utils.formats_data import format_temporal_variable,define_output_path



def operadar(filename:str,
             modelname:str=cf.model,
             read_tmatrix:bool=True,
             in_dir_path:str=cf.input_filePath,
             out_dir_path:str=cf.outPath,
             tmatrix_path:str=cf.path_Tmatrix,
             microphysics_scheme:str=cf.micro_scheme,
             hydrometeorMoments:dict[int]=cf.hydrometeors_moments,
             radar_band:str=cf.radar_band,
             radarloc:str|list=cf.radarloc,
             distmax_rad:float=cf.distmax_rad,
             Tmatrix_params:dict={},
             mixed_phase_parametrization:str=cf.MixedPhase,
             subDomain:list[float]|None=cf.subDomain,
             get_more_details=False,
             append_in_file=False,
             ) -> tuple[bool,bool,dict]:
    """Radar forward operator main function. 

    Args:
        filename (str): Only the name of the file.
        modelname (str): Defaults to cf.model.
        read_tmatrix (bool, optional): Option to save computing time within a loop.
                Need to be set to True at least for the first iteration.
                Defaults to True.
        in_dir_path (str, optional): Path where the input files are stored.
                Defaults to cf.input_filePath.
        out_dir_path (str, optional): Path to store the output files.
                Defaults to cf.outPath.
        tmatrix_path (str, optional) : Path where the Tmatrix lookup tables are stored.
        microphysics_scheme (str, optional): name of the microphysics scheme.
                Defaults to cf.micro_scheme.
        hydrometeorMoments (dict of form {str : int}), optional): dict containing the number of moments
                for each hydrometeor of the microphysics scheme.
                Defaults to cf.hydrometeors_moments
        radar_band (str, optional): Defaults to cf.radar_band.
        radarloc (str | list): location of the radar for simulation.
                Either 'center' (i.e. center of the grid) or a [lat_radar,lon_radar]
                coordinate. Defaults to cf.radarloc.
        distmax_rad (float): Maximum radius of the radar data to compute pseudo-observations.
                DEfaults to cf.distmax_rad.
        Tmatrix_params (dict, optional): dictionary containing the Tmatrix tables parameters.
                Defaults to {}.
        mixed_phase_parametrization (str, optional): Defaults to cf.MixedPhase.
        subDomain (list[float] | None, optional): Defaults to cf.subDomain.
        get_more_details (bool): Defaults to False.
        append_in_file (bool): Defaults to False.

    Returns:
        read_tmatrix (bool): to update the boolean during a loop
        Tmatrix_params (dict) : to keep in memory the Tmatrix parameters and values throughout multiple iterations over the same radar band
        
    """
    
    begin_program_timer = tm.time()
    
    if append_in_file and subDomain != None :
        print('Cannot reinject a subdomain into the original input file. Please, set subDomaine = None to continue.')
        exit()
    print(cf.real_case)
    # Create or check tree structure of the output directory path
    create_tree_structure_outFiles(output_dir=Path(out_dir_path))
    # Format temporal variable and output file name
<<<<<<< HEAD
    input_file_path = Path(in_dir_path+filename)    
=======
    input_file_path = Path(in_dir_path+filename)
>>>>>>> 7359c82 (Integration of BenoitVie modifications)
    temporal_variable = format_temporal_variable(filePath=input_file_path,
                                                 model_type=modelname,
						 real_case=cf.real_case,
                                                 )
    outFilePath = define_output_path(out_dir_path=out_dir_path,
                                     model=modelname,
                                     scheme=microphysics_scheme,
                                     radar_band=radar_band,
                                     temporal_variable=temporal_variable,
                                     ) 
    if not Path(outFilePath).with_suffix('.nc').exists() or append_in_file :
        
        # Read Tmatrix tables (files from Clotilde)
        if read_tmatrix :
            Tmatrix_hydromet_list = link_keys_with_available_hydrometeors(hydrometeorMoments=hydrometeorMoments,
                                                                          datatype='tmatrix',
                                                                          )
            Tmatrix_params = read_Tmatrix_Clotilde(band=radar_band,
                                                   scheme=microphysics_scheme,
                                                   pathTmat=tmatrix_path,
                                                   hydrometeors=Tmatrix_hydromet_list,
                                                   verbose=get_more_details,
                                                   )
        # Read model variables
        [X, Y, Alt, lon, lat, M, Nc, Tc] = read_model_file(filePath=input_file_path,
                                                           modelname=modelname,
                                                           domain=subDomain,
                                                           hydrometeorMoments=hydrometeorMoments,
                                                           verbose=get_more_details,
                                                           )
        # Compute radar geometry
        mask_dist_max, elevations = compute_radar_geometry(X=X, Y=Y, Z=Alt, Tc=Tc,
                                                           elev_max=Tmatrix_params['ELEVmax']["rr"],
                                                           model=modelname,
                                                           distmax_rad=distmax_rad,
                                                           radarloc=radarloc,
                                                           )
        # Mask precipitations
        mask_precip = mask_precipitations(contents=M,
                                          expMmin=Tmatrix_params['expMmin']["rr"],
                                          hydrometeors_moments=hydrometeorMoments,
                                          )
        # Combine masks
        partial_mask = (mask_precip & mask_dist_max)
        
        # Compute mixed phase parametrization
        [M, Nc,Fw] = compute_mixed_phase(contents=M,
                                         concentrations=Nc,
                                         hydrometeorMoments=hydrometeorMoments,
                                         expMmin=Tmatrix_params['expMmin']["rr"],
                                         parametrization=mixed_phase_parametrization,
                                         ) 
        # Compute dual-pol radar variables
        dpolDict = compute_dualpol_variables(temperature=Tc,
                                             mask_precip_dist=partial_mask,
                                             elev=elevations, Fw=Fw,
                                             contents=M,
                                             concentrations=Nc,
                                             tmatrix_param=Tmatrix_params,
                                             hydrometeorMoments=hydrometeorMoments,
                                             X=X, Y=Y, Z=Alt,
                                             lat=lat, lon=lon,
                                             date_time=temporal_variable,
                                             output_file_path=outFilePath,
                                             append_in_fa=append_in_file,
                                             )
        # Saving file or reinjecting fields into the input file
        if append_in_file :
            del M, Nc, Fw, Alt, lat, lon, Tc, elevations
            append_in_input_file(complete_input_path=input_file_path,
                                 dpolVar=dpolDict,
                                 var2add=cf.dpol2add,
                                 )
        else :
            save_netcdf(X=X, Y=Y, Z=Alt, lat=lat, lon=lon,
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
        read_tmatrix = True
        print("File exists at :",outFilePath+'.nc')
        print("-----------------------------------------------------------------------------")
        return read_tmatrix, dict()




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="OPERADAR")
    parser.add_argument("filename", type=str)
    parser.add_argument("--append", action="store_true",default=False)
    parser.add_argument("--verbose", action="store_true",default=False)
    operad_args = parser.parse_args()
    operadar(filename=operad_args.filename,
             get_more_details=operad_args.verbose,
             append_in_file=operad_args.append)
