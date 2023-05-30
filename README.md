# operadar
Computes dual-pol variables (Zh, Zdr, Kdp, Rhohv) in the 3D model grid for Arome or MesoNH model using existing Tmatrix tables
input: Tmatrix tables and model file (Arome fa or MesoNH netcdf)
output: netcdf file with lat,lon (or X, Y) + Zh, Zdr, Kdp, Rhohv, T, Alt 

1) create a configuration file like operad_conf_AROME_ICE3.py or operad_conf_MesoNH_ICEidpc.py
==> select your model, microphysics scheme options
!!! ICE3 => list_types=['vv','cc','rr','ii','ss','gg']
            list_types_tot=['rr','ii','ss','gg','wg']
!!! ICE4 => list_types=['vv','cc','rr','ii','ss','gg','hh']
            list_types_tot=['rr','ii','ss','gg','wg','wh']

MesoNH only:
LIMA => list_types=['vv','cc','rr','ii','ss','gg']
            list_types_tot=['rr','ii','ss','gg','wg']
LIMH => list_types=['vv','cc','rr','ii','ss','gg','hh']
            list_types_tot=['rr','ii','ss','gg','wg','wh']

==> select radar option (band)

==> change list of time steps, Tmatrix table and model file directories

2) change the name of the config file in operad.py

3) python3 operad.py

===================================
@todo
- find a way to read the ice concentration in AROME file ?
- implement the LIMA/LIMH options for AROME (need rain and cloud water concentration) 
- add the radar geometry option (with elevations and beam filtering with gaussian)
- link this project with tmatrix DPOLSIMUL git (which produces the required tables !)
