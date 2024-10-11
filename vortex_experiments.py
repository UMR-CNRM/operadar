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