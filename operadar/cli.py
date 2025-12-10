#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from .forward_operator import operadar
from .load_config_file import load_configuration_file


def run_cli():
   	# Reading input arguments
    parser = argparse.ArgumentParser(description="Observation oPErator for polarimetric RADAR (OPERADAR). This software is developped at the CNRM, France.")
    parser.add_argument("filename", type=str, help="Only the name of the input file, not the entire path. Please use the config file to provide the path to access the file.")
    parser.add_argument("config",   type=str, help="Name of the configuration file. Before running the code, you need to create a configuration file in ./configFiles/ based on the template provided.")
    parser.add_argument("--append", action="store_true", default=False, help="Append the computed dual-polarimetric fields into the input file (default value: False).")
    parser.add_argument("--verbose", action="store_true", default=False, help="If argument provided, show more details if activated (default value: False).")
    operad_args = parser.parse_args()
    
    config = load_configuration_file(operad_args.config)
    operadar(filename=operad_args.filename,
             configuration=config,
             get_more_details=operad_args.verbose,
             append_in_file=operad_args.append)
