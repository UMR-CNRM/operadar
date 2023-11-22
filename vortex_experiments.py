def set_vortex_experiments(run: str, microphysics_scheme: str) :
    expeOLIVE = ''
    if microphysics_scheme == 'ICE3' :
        if run == '00' :
            expeOLIVE = 'GN49'
        elif run == '12' :
            expeOLIVE = 'GOIM'
    elif microphysics_scheme == 'ICE4' :
        if run == '00' :
            expeOLIVE = 'GN51'
        elif run == '12' :
            expeOLIVE = 'GOJI'
    elif microphysics_scheme == 'LIMASG' :
        if run == '00' :
            expeOLIVE = 'GOTW'
        elif run == '12' :
            expeOLIVE = 'GOVA'
    elif microphysics_scheme == 'LIMAAG' :
        if run == '00' :
            expeOLIVE = 'GOIP'
        elif run == '12' :
            expeOLIVE = 'GOVE'
    return expeOLIVE