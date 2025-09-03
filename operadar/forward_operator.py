#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 13:17:47 2018
@author: augrosc & lebastardt & montangonm & sinhorin & davidcl
"""

# External modules
import os
import time as tm
from sys import exit
from pathlib import Path
from dataclasses import dataclass
from typing import List, Union, Optional, Dict

# 0PERADAR modules
from operadar.read.model import read_model_file
from operadar.utils.masking import mask_precipitations
from operadar.radar.geometry import compute_radar_geometry
from operadar.save.append_dpolvar import append_in_input_file
from operadar.microphysics.mixed_phase import compute_mixed_phase
from operadar.radar.dualpol_variables import compute_dualpol_variables
from operadar.read.lookup_tables import read_and_extract_tables_content
from operadar.utils.make_links import link_keys_with_available_hydrometeors
from operadar.save.save_dpolvar import create_tree_structure_outFiles, save_netcdf
from operadar.utils.formats_data import format_temporal_variable, define_output_path



@dataclass
class Config:
    input_filePath: str
    output_filePath: str
    path_tables: str
    model: str
    real_case: bool
    microphysics_scheme: str
    hydrometeors_moments: Dict[str, int]
    cloud_water_over: str
    subDomain: Optional[List[Union[float,int]]]
    mixed_phase_parametrization: str
    save_netcdf_single_hydrometeor: bool
    dpol2add: List[str]
    scattering_method: str
    radar_band: str
    distmax_rad: float
    radarloc: Optional[str]
    cnst_angle: float



def load_config_instance(config_file) -> Config:
    print('user file config',config_file)
    return Config(**{key: val for key,val in vars(config_file).items() if not key.startswith("__") and not callable(val)})



def update_config_with_overrides(config, **overrides):
    for key, value in overrides.items():
        if value is not None and hasattr(config, key):
            setattr(config, key, value)



def operadar(filename: str,
             configuration:str,
             read_tables: bool=True,
             in_dir_path=None,
             out_dir_path=None,
             tables_path=None,
             modelname=None,
             real_case=None,
             microphysics_scheme=None,
             hydrometeorMoments=None,
             cloud_water_over=None,
             subDomain=None,
             mixed_phase_parametrization=None,
             save_netcdf_single_hydrometeor=None,
             dpol2add=None,
             scattering_method=None,
             radar_band=None,
             distmax_rad=None,
             radarloc=None,
             cnst_angle=None,
             tables_content:dict={},
             get_more_details=False,
             append_in_file=False,
             ) -> tuple[bool,dict]:
    """Radar forward operator main function. 

    Args:
        filename (str): Only the name of the file.
        modelname (str): Defaults to cf.model.
        read_tables (bool, optional): Option to save computing time within a loop.
                Need to be set to True at least for the first iteration.
                Defaults to True.
        in_dir_path (str, optional): Path where the input files are stored.
                Defaults to cf.input_filePath.
        out_dir_path (str, optional): Path to store the output files.
                Defaults to cf.outPath.
        tables_path (str, optional) : Path where the lookup tables are stored.
        microphysics_scheme (str, optional): name of the microphysics scheme.
                Defaults to cf.micro_scheme.
        hydrometeorMoments (dict of form {str : int}), optional): dict containing the number of moments
                for each hydrometeor of the microphysics scheme.
                Defaults to cf.hydrometeors_moments
        radar_band (str, optional): Defaults to cf.radar_band.
        radarloc (str | list | None) : location of the radar for simulation.
                Either 'center' (i.e. center of the grid) or a [lat_radar,lon_radar]
                coordinate or None. Defaults to cf.radarloc.
        distmax_rad (float): Maximum radius of the radar data to compute pseudo-observations.
                DEfaults to cf.distmax_rad.
        tables_content (dict, optional): dictionary containing the tables parameters and required columns.
                Defaults to {}.
        mixed_phase_parametrization (str, optional): Defaults to cf.MixedPhase.
        subDomain (list[float] | None, optional): Defaults to cf.subDomain.
        get_more_details (bool): Defaults to False.
        append_in_file (bool): Defaults to False.

    Returns:
        read_tables (bool): to update the boolean during a loop
        tables_content (dict) : to keep in memory the tables parameters and values throughout multiple iterations over the same radar band
        
    """
    begin_program_timer = tm.time()
    
    #Load the configuration based on the config file or the given arguments
    conf = load_config_instance(config_file=configuration)
    update_config_with_overrides(conf,
                                 input_filePath=in_dir_path, output_filePath=out_dir_path, path_tables=tables_path,
                                 model=modelname, real_case=real_case, microphysics_scheme=microphysics_scheme,
                                 hydrometeors_moments=hydrometeorMoments, cloud_water_over=cloud_water_over,
                                 subDomain=subDomain, mixed_phase_parametrization=mixed_phase_parametrization,
                                 save_netcdf_single_hydrometeor=save_netcdf_single_hydrometeor, 
                                 dpol2add=dpol2add, scattering_method=scattering_method, radar_band=radar_band,
                                 distmax_rad=distmax_rad, radarloc=radarloc, cnst_angle=cnst_angle)
    
    
    if append_in_file and subDomain != None :
        print('Cannot reinject a subdomain into the original input file. Please, set subDomaine = None to continue.')
        exit()
    
    # Create or check tree structure of the output directory path
    create_tree_structure_outFiles(output_dir=Path(conf.output_filePath))
    # Format temporal variable and output file name
    input_file_path = Path(os.path.join(conf.input_filePath,filename))
    temporal_variable = format_temporal_variable(filePath=input_file_path,
                                                 model_type=conf.model,
						                         real_case=conf.real_case,
                                                 )
    outFilePath = define_output_path(out_dir_path=conf.output_filePath,
                                     model=conf.model,
                                     scheme=conf.microphysics_scheme,
                                     radar_band=conf.radar_band,
                                     temporal_variable=temporal_variable,
                                     ) 
    if not outFilePath.with_suffix('.nc').exists() or append_in_file :
        
        # Read lookup tables
        if read_tables :
            hydromet_list = link_keys_with_available_hydrometeors(hydrometeorMoments=conf.hydrometeors_moments,
                                                                  datatype='tables',
                                                                  )
            tables_content = read_and_extract_tables_content(band=conf.radar_band,
                                                             scheme=conf.microphysics_scheme,
                                                             moments=conf.hydrometeors_moments,
                                                             path_table=conf.path_tables,
                                                             dpol2add=conf.dpol2add,
                                                             hydrometeors=hydromet_list,
                                                             cloud_water_over=conf.cloud_water_over,
                                                             verbose=get_more_details,
                                                            )
        # Read model variables
        [X, Y, Alt, lon, lat, M, Nc, Tc] = read_model_file(filePath=input_file_path,
                                                           modelname=conf.model,
                                                           micro_scheme=conf.microphysics_scheme,
                                                           real_case=conf.real_case,
                                                           domain=conf.subDomain,
                                                           hydrometeorMoments=conf.hydrometeors_moments,
                                                           verbose=get_more_details,
                                                           )
        # Compute radar geometry
        mask_dist_max, elevations = compute_radar_geometry(X=X, Y=Y, Z=Alt, Tc=Tc,
                                                           elev_max=tables_content['ELEVmax']["rr"],
                                                           model=conf.model,
                                                           cnst_angle=conf.cnst_angle,
                                                           distmax_rad=conf.distmax_rad,
                                                           radarloc=conf.radarloc,
                                                           )
        # Mask precipitations
        mask_precip = mask_precipitations(contents=M,
                                          expMmin=tables_content['expMmin']["rr"],
                                          hydrometeors_moments=conf.hydrometeors_moments,
                                          )
        # Combine masks
        partial_mask = (mask_precip & mask_dist_max)
        
        # Compute mixed phase parametrization
        [M, Nc,Fw] = compute_mixed_phase(contents=M,
                                         concentrations=Nc,
                                         hydrometeorMoments=conf.hydrometeors_moments,
                                         expMmin=tables_content['expMmin']["rr"],
                                         parametrization=conf.mixed_phase_parametrization,
                                         ) 
        # Compute dual-pol radar variables
        dpolDict = compute_dualpol_variables(temperature=Tc,
                                             mask_precip_dist=partial_mask,
                                             elev=elevations,
                                             Fw=Fw,
                                             contents=M,
                                             concentrations=Nc,
                                             tables_dict=tables_content,
                                             X=X, Y=Y, Z=Alt,
                                             lat=lat, lon=lon,
                                             date_time=temporal_variable,
                                             output_file_path=outFilePath,
                                             append_in_fa=append_in_file,
                                             config=conf,
                                             )
        
        # Saving file or reinjecting fields into the input file
        if append_in_file :
            del M, Nc, Fw, Alt, lat, lon, Tc, elevations
            append_in_input_file(complete_input_path=input_file_path,
                                 dpolVar=dpolDict,
                                 var2add=conf.dpol2add,
                                 )
        else :
            save_netcdf(X=X, Y=Y, Z=Alt, lat=lat, lon=lon, datetime=temporal_variable,
                        dpolDict=dpolDict, contentsDict=M, concentrationsDict=Nc,
                        temperature=Tc, outfile=outFilePath, config=conf,
                        )
        # For multiple iterations over different time but with the same settings, save time by
        # not reading again lookup tables and lat/lon fields (if available)
        read_tables = False
        elapsed_time = tm.time() - begin_program_timer
        print("Elapsed time :",int(elapsed_time//60),"minutes",int(elapsed_time%60),"seconds")
        print("-----------------------------------------------------------------------------")
        return read_tables, tables_content
    
    else :
        read_tables = True
        print("File exists at :",outFilePath.with_suffix('.nc'))
        print("-----------------------------------------------------------------------------")
        return read_tables, dict()

