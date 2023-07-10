# operadar
Computes dual-pol variables (Zh, Zdr, Kdp, Rhohv) in the 3D model grid for Arome or MesoNH model using existing Tmatrix tables
* INPUT  : Tmatrix tables and model file (Arome fa or MesoNH netcdf)
* OUTPUT : netcdf file with lat,lon (or X, Y) + Zh, Zdr, Kdp, Rhohv, T, Alt 

## How to run operadar
1) Create a configuration file (like `operad_conf_AROME_ICE3.py` or `operad_conf_MesoNH_ICEidpc.py`) for the chosen model
2) Select your model microphysics scheme options
   - ICE3 or LIMA (without hail) :
     - `htypes_model   = ['vv','cc','rr','ii','ss','gg']` # model's hydrometeors related variables
     - `list_types_tot = ['rr','ii','ss','gg','wg']`    # idem **+ wet hydrometeors (computed)**
   - ICE4 or LIMA with hail :
     - htypes_model=['vv','cc','rr','ii','ss','gg','hh']
     - list_types_tot=['rr','ii','ss','gg','wg','hh','wh']
3) Select a radar band (C, S or X)
4) Change time steps, Tmatrix table and model file directories/paths
5) Run in a terminal `python3 operad.py`


## To Do
- [ ] find a way to read the ice concentration in AROME file ?
- [ ] implement the LIMA/LIMH options for AROME (need rain and cloud water concentration) 
- [ ] add the radar geometry option (with elevations and beam filtering with gaussian)
- [ ] link this project with tmatrix DPOLSIMUL git (which produces the required tables !)
