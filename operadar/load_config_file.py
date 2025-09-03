#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import importlib.util as iu

    
OPERADAR_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIGS_DIR = os.path.join(OPERADAR_DIR,"../configFiles")    
    
def load_configuration_file(user_file_name: str):
    if user_file_name[-3:]==".py":
    	filename = user_file_name
    else :
        filename = f"{user_file_name}.py"
    path = os.path.join(CONFIGS_DIR, filename)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Cannot find : {path}")    
    spec = iu.spec_from_file_location(user_file_name, path)
    config_module = iu.module_from_spec(spec)
    spec.loader.exec_module(config_module)
    return config_module
