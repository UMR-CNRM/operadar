#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def get_vortex_experiments(csvRow: str, microphysics_scheme: str) :
    listExpe = [str(x) for x in csvRow.expeNames.strip().split(',')]
    if microphysics_scheme == 'ICE3' :
        return listExpe[0]
    elif microphysics_scheme == 'ICE4' :
        return listExpe[1]
    elif microphysics_scheme == 'LIMASG' :
        return listExpe[2]
    elif microphysics_scheme == 'LIMAAG' :
        return listExpe[3]
    elif microphysics_scheme == 'LIMA49t' :
        return listExpe[4]