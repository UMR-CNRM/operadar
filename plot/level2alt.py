import numpy as np
import xarray as xr
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter


def get_iso0_alt(ds):
    ds['iso0'] = ds.level.where((ds.T <= 1)&(ds.T >= -1)).idxmin('level')
    ds['alt_iso0']=ds.Alt.where(ds.iso0==ds.level,0)
    return ds.alt_iso0.sum('level')


def interpolate_variables(ds,resolV=500,alti_min=0,alti_max=15e3):
    print("  Interpolation from model levels to meter altitude")
    # Extract values and the irregular coordinates (x,y,z) then flatten
    x = ds.x.expand_dims({"level":len(ds.level),"y":len(ds.y)}).values.ravel()
    y = ds.y.expand_dims({"level":len(ds.level),"x":len(ds.x)},axis=(0,2)).values.ravel()
    z = ds.Alt.values.ravel()
    zh = ds.Zh.values.ravel()
    kdp = ds.Kdp.values.ravel()
    zdr = ds.Zdr.values.ravel()
    rhoHV = ds.Rhohv.values.ravel()
    # ccRain = ds.CCrain.values.ravel()
    # ccIce = ds.CCice.values.ravel()

    # Creating grid with regular x,y and z 
    Y, Z, X = np.meshgrid(ds.y.values,np.arange(alti_min,alti_max+resolV,resolV),ds.x.values)

    # Interpolation of the dual-pol variables
    zh_interp  = griddata((z,y,x), zh, (Z,Y,X), method='nearest')
    zdr_interp = griddata((z,y,x), zdr, (Z,Y,X), method='nearest')
    kdp_interp = griddata((z,y,x), kdp, (Z,Y,X), method='nearest')
    rhoHV_interp = griddata((z,y,x), rhoHV, (Z,Y,X), method='nearest')
    
    # # Interpolation of hydrometeors concentrations
    # ccRain_interp = griddata((z,y,x), ccRain, (Z,Y,X), method='nearest')
    # ccIce_interp = griddata((z,y,x), ccIce, (Z,Y,X), method='nearest')
    
    # # Interpolation of hydrometeor contents
    # hydro_dim = len(ds.hydrometeor.values)
    # M_interp = np.zeros((hydro_dim,Y.shape[0],Y.shape[1],Y.shape[2]))
    # for i in range(hydro_dim) :
    #     extract_hydrometeor_data = ds.M[i].values.ravel()
    #     interp_hydrometeor_data = griddata((z,y,x), extract_hydrometeor_data , (Z,Y,X), method='nearest')
    #     M_interp[i] += interp_hydrometeor_data
    
    return zh_interp, zdr_interp, kdp_interp, rhoHV_interp #, ccRain_interp, ccIce_interp, M_interp


def interpolate_dataset(ds,lvl_intervals=1,resolV=500,alti_min=0,alti_max=15e3):
    # Reconstruction du champ 2D d'altitude de l'iso0
    #ds['alt_iso0']=get_iso0_alt(ds)
    
    sub_ds = ds.sel(level=slice(9,89,lvl_intervals))
    sub_ds.coords['x'] = sub_ds.x.astype('int32')*1000 
    sub_ds.coords['y'] = sub_ds.y.astype('int32')*1000

    #zh_interp, zdr_interp, kdp_interp, rhoHV_interp, ccRain_interp, ccIce_interp, M_interp = interpolate_variables(sub_ds,cfg_dict)
    zh_interp, zdr_interp, kdp_interp, rhoHV_interp = interpolate_variables(sub_ds,resolV,alti_min,alti_max)
    
    # Storage in a new interpolated dataset
    ds_model = xr.Dataset(

    data_vars=dict(
        zh=(["z","y","x"], zh_interp.astype("float32"),{"units":"dBZ"}),
        zdr=(["z","y","x"], zdr_interp.astype("float32"),{"units":"dB"}),
        kdp=(["z","y","x"], kdp_interp.astype("float32"),{"units":"°/km"}),
        rhohv=(["z","y","x"], rhoHV_interp.astype("float32"),{"units":"1"}),
        # cc_ice=(["z","y","x"],ccIce_interp.astype("float32"),{"units":"nb/m3"}),
        # cc_rain=(["z","y","x"],ccRain_interp.astype("float32"),{"units":"nb/m3"}),
        # M=(["hydrometeor","z","y","x"],M_interp.astype("float32"),{"units":"g/dry air kg"}),
        #iso0=(["y","x"], gaussian_filter(sub_ds.alt_iso0.values ,sigma=10).astype("float32"),{"units":"m"}),
    ),
    coords=dict(
        x =(["x"],sub_ds.x.values.astype("int32"),{"units":"m"}),
        y =(["y"],sub_ds.y.values.astype("int32"),{"units":"m"}),
        lon=(["y", "x"], sub_ds.lon.values.astype("float32")),
        lat=(["y", "x"], sub_ds.lat.values.astype("float32")),
        z=(["z"],np.arange(alti_min,alti_max+resolV,resolV).astype("int16"),{"units":"m"}),
        #hydrometeor = (["hydrometeor"],sub_ds.hydrometeor.values),
        time=(ds.time.values.astype("datetime64"))
    ),
    )
    
    return ds_model

