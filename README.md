# operadar
Computes dual-pol variables (Zh, Zdr, Kdp, Rhohv) in the 3D model grid for Arome or MesoNH model using existing Tmatrix tables
INPUT  : Tmatrix tables and model file (Arome fa or MesoNH netcdf)
OUTPUT : netcdf file with lat,lon (or X, Y) + Zh, Zdr, Kdp, Rhohv, T, Alt 

## How to run operadar
1) Create a configuration file (like `operad_conf_AROME_ICE3.py` or `operad_conf_MesoNH_ICEidpc.py`) for the chosen model
2) Select your model microphysics scheme options
   - <p>ICE3 => htypes_model=['vv','cc','rr','ii','ss','gg'] # available model variables<br>
           list_types_tot=['rr','ii','ss','gg','wg']</p>
   - ICE4 => htypes_model=['vv','cc','rr','ii','ss','gg','hh']
            list_types_tot=['rr','ii','ss','gg','wg','hh','wh']

LIMA => htypes_model=['vv','cc','rr','ii','ss','gg']
            list_types_tot=['rr','ii','ss','gg','wg']
LIMH => htypes_model=['vv','cc','rr','ii','ss','gg','hh']
            list_types_tot=['rr','ii','ss','gg','wg','hh','wh']

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
