# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 13:20:49 2015

@author: chwala-c
"""

import numpy as np
#from pkg_resources import resource_filename
from importlib.resources import files

def MPM(f, P, T, U, wa, wae, R, output_type='ref'):
    """
    Calculate complex refractivity of the atmosphere
    
    Use the MPM (millimeter wave propagation model) to calculate the 
    refractivity of the atmosphere at a given frequency. The output 
    can be automatically transformed from refractivity to more usabel
    quantities.
    
    Note that `wa`, `wae` and `R` are not yet supported. They can be provided,
    but do not affect the calculation.
    
    Parameters
    ----------
    
    f : int, float, or list or array of these
        Frequency in GHz
    P : float
        Air pressure in mbar
    T : float
        Air temperature in degree Celcius
    U : float
        Relative humidity in %
    wa : float
        Density of water droplets in g/m^3
    wae : float
        Concentration of ice particles in g/m^3
    R : float
        Rain rate
    output_tpye : str, optional
        Desired conversion of complex refractivity for output. Default is
        `ref` which returns the complex refractiviy
        Supported types are:
        'ref' = Refractivity
        'att' = Attenuation in db/km
        'dis' = Phase dispersion in deg/km
        'del' = Group delay in ps/km
        'abs' = Absorption coefficients in 1/m
        
    Returns
    -------
    
    float, np.array
        The desired conversion of the compelx refractiviy
    
    """
    print('MPM')
    
    param = inputconv(P, T, U)
    
    NV = watervapormodule(f, param.e, param.pd, param.th)
    ND = dryairmodule(f, param.e, param.pd, param.th)
    
    # Sum up all refractivities
    N = NV + ND
    
    return outconv(f, N, output_type)
    
def inputconv(P, T, U):
    """
    Convert input parameters P, T, U to e, pd, th

    Parameters
    ----------
    P : float or numpy.ndarray
        Air pressure in mbar
    T : float or numpy.ndarray
        Air temperature in degree Celcius
    U : float or numpy.ndarray
        Relative humidity in %

    Returns
    -------
    a : namedtuple
        Namedtuple of
         * e (Partial water vapor pressure)
         * pd (Partial pressure of dry air)
         * th (Reciprocal temperature)
    """

    from collections import namedtuple
    a = namedtuple('a', ['e', 'pd', 'th'])

    # Convert to numpy array if they aren't already
    P = np.asarray(P)
    T = np.asarray(T)
    U = np.asarray(U)

    # Check for invalid temperatures
    if np.any((T + 273.15) <= 0):
        raise ValueError("Absolute temperature must be positive (T > -273.15°C)")

    TK = T + 273.15
    # Reciprocal temperature theta
    th = 300.0/TK
    # Calculate log of water vapor saturation pressure
    log_es = np.log(2.408e11) + 5*np.log(th) - 22.644*th

    # Handle potential underflow in exp() calculation
    es = np.zeros_like(log_es)
    underflow_mask = log_es < -700  # exp(-700) ~ 1e-304, near underflow
    es[~underflow_mask] = np.exp(log_es[~underflow_mask])

    # Partial pressure of water vapor
    e = es * U / 100.0
    # Partial pressure of dry air
    pd = P - e
    # water vapor density in g/m3
    Rs = 461.25  # Specific gas constant for water vapor (Stöcker, Taschenbuch der Physik)
    rhogas = 100000 * e / (Rs * TK)

    return a(e=e, pd=pd, th=th)
    
def outconv(f, N, output_type):
    """ 
    Convert complex refractivity N to desired output 
    
    Parameters
    ----------
    
    f : float, np.array of float
        Frequency in GHz
    N : complex, np.array of complex
        Refractivity
    output_type: str
        Desired output tpye. Supported types are:
        'ref' = Refractivity
        'att' = Attenuation in db/km
        'dis' = Phase dispersion in deg/km
        'del' = Group delay in ps/km
        'abs' = Absorption coefficients in 1/m
    
    Returns
    -------
    
    out : float, np.array
        The desired conversion of the refractiviy
    """
    
    c=299792458
    
    if len(N) == 0:
        return 0
    else:
        if output_type == 'ref':
            # Refractivity
            out = N
        elif output_type == 'att':
            # Attenuation in db/km
            out = 0.1820*f*np.imag(N)
        elif output_type == 'dis':
            # Phase dispersion beta in deg/km
            out = 1.2008*f*np.real(N)
        elif output_type == 'del':
            # Group delay tay in ps/km
            out = 3.3356 * np.real(N)
        elif output_type == 'abs':
            # Absorption coefficients in 1/m
            out = 4*np.pi*1000/c*f*np.imag(N)
        else:
            raise ValueError('output_type not supported')
        return out
          
def dryairmodule(f_vec, e, pd, th):
    """
    Calculate refractivity for oxygen lines
    
    Parameters
    ----------
    
    blabla
    
    """
    
    # # Make sure f_vec is an iterabel array
    # if type(f_vec) == float or type(f_vec) == int:
    #     f_vec = np.array([f_vec,])
    # Make sure f_vec is an iterable array (handles both standard types and numpy.float64)
    f_vec = np.atleast_1d(f_vec)
    
    ND = []
    # Note: We use the full package path 'operadar.radar.pyMPM', not just 'pyMPM'
    package_path = files('operadar.radar.pyMPM')
    with package_path.joinpath('oxygen93.txt').open('r') as f:
        oxygen_lines = np.loadtxt(f)
        #oxygen_lines = np.loadtxt(resource_filename('pyMPM',
        #                                        'oxygen93.txt'))
    
    for f in f_vec:
        # Non dispersive part
        nd = 0.2588*pd*th
        # Loop over lines
        for line in oxygen_lines:
            f_line = line[0]
            # Line strength (with correction factor 1e-6)
            S =  1e-6*line[1]/f_line*pd*th**3*np.exp(line[2]*(1-th))
            # Line width
            gamma = line[3]/1000.0*(pd*th**line[4] + 1.1*e*th)
            # Zeeman effect
            gamma = np.sqrt(gamma**2 + 2.25e-6)
            # Overlap parameter (with correction factor 1e-3)
            delta = 1e-3*(line[5] + line[6]*th)*(pd+e)*th**0.8
            # Line form function
            F = f*((1-1j*delta)/(f_line-f-1j*gamma) - (1+1j*delta)/(f_line+f+1j*gamma))
            nd = nd + S*F
        # Non resonant part
        So = 6.14e-5*pd*th**2
        Fo = -f/(f+1j*0.56e-3*(pd+e)*th**0.8)
        Sn = 1.4e-12*pd**2*th**3.5
        Fn = f/(1 + 1.93e-5*f**1.5) # Correction 1.9 -> 1.93
        nd = nd + So*Fo + 1j*Sn*Fn
        ND.append(nd)
    return np.array(ND)

        
def watervapormodule(f_vec, e, pd, th):
    """ 
    Calcualte refractivity for water vapor lines 
    
    Parameters
    ----------
    
    f_vec : float, or array of floats
        Frequency in GHz
        
        ....
        
    """
    
    # # Make sure f_vec is an iterabel array
    # if type(f_vec) == float or type(f_vec) == int:
    #     f_vec = np.array([f_vec,])
    # Make sure f_vec is an iterable array (handles both standard types and numpy.float64)
    f_vec = np.atleast_1d(f_vec)
    
    NV = []
    #water_lines = np.loadtxt(resource_filename('pyMPM',
    #                                           'water93.txt'))
    # Note: We use the full package path 'operadar.radar.pyMPM', not just 'pyMPM'
    package_path = files('operadar.radar.pyMPM')
    with package_path.joinpath('water93.txt').open('r') as f:
        water_lines = np.loadtxt(f)
    
    # Loop over frequencies
    for f in f_vec:
        nv = (4.163*th + 0.239)*e*th
        # Loop over lines
        for line in water_lines:
            f_line = line[0]
            # Line strength
            S = line[1]/f_line*e*th**(3.5)*np.exp(line[2]*(1 - th))
            # Line width
            gamma = line[3]/1000 * (line[4]*e*th**line[6] + pd*th**line[5])
            # Doppler broadening
            #i_dopp = np.where(pd + e < 0.7)
            # avant : if pd + e < 0.7:

            ind = np.where((pd + e )<0.7)
            gamma[ind] = (0.535*gamma + np.sqrt(0.217*gamma**2 + (1.46e-6*gamma*np.sqrt(th))**2))[ind]
            # Line form function
            F = f*(1/(f_line - f - 1j*gamma) - 1/(f_line + f + 1j*gamma))
            nv = nv + S*F
        NV.append(nv)
    return np.array(NV)

# ======================================================================


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # -------------------------------
    # Demo: Ground radar vertical cross-section
    # -------------------------------

    # Grid dimensions: (height_levels, range_gates)
    Nz = 50      # height levels (up to ~10 km)
    Nx = 100     # range gates (up to 50 km)
    Ny = 1       # single azimuth (2D cross-section); we'll squeeze later

    # Height grid (meters) – exponential spacing or linear
    z = np.linspace(0, 10000, Nz)  # 0 to 10 km

    # Range grid (meters)
    r = np.linspace(0, 50000, Nx)  # 0 to 50 km

    # Create 2D mesh (Nz, Nx), then expand to (Nz, 1, Nx) for 3D compatibility
    Z, R = np.meshgrid(z, r, indexing='ij')  # Z[z_idx, r_idx], R[z_idx, r_idx]

    # --- Build atmospheric profiles ---
    # Standard atmosphere: temperature decreases with height (~6.5 K/km)
    T_K = 288.15 - 0.0065 * Z  # K
    T_C = T_K - 273.15          # °C

    # Pressure: barometric formula (approx)
    P_Pa = 101325 * (T_K / 288.15) ** (9.80665 * 0.0289644 / (8.31447 * 0.0065))
    P_mbar = P_Pa / 100.0       # Convert Pa → mbar (1 mbar = 100 Pa)

    # Relative humidity: decreases with height (e.g., 80% at surface, 20% at 5 km)
    RH = 80.0 * np.exp(-Z / 2000.0)  # %
    RH = np.clip(RH, 1.0, 100.0)     # Avoid <1%

    # Expand to 3D: (Nz, Ny=1, Nx)
    pressure_3d = P_mbar[:, np.newaxis, :]      # (Nz, 1, Nx)
    temperature_3d = T_K[:, np.newaxis, :]      # in K → will convert to °C in call
    rh_3d = RH[:, np.newaxis, :]

    # Radar frequency (GHz)
    freq_ghz = 94.0  # W-band cloud radar

    # Compute attenuation as in your usage
    try:
        kext_gaz_raw = (
            1e-3 / 4.343
            * MPM(
                freq_ghz,
                pressure_3d * 1e-2,           # Pa → hPa (since P_mbar = Pa/100, this gives hPa = mbar)
                temperature_3d - 273.15,      # K → °C
                rh_3d,
                wa='None',
                wae='None',
                R='None',
                output_type='att'             # dB/km
            )
        )
    except Exception as e:
        print(f"⚠️ MPM failed: {e}")
        print("Make sure 'oxygen93.txt' and 'water93.txt' are in the package data directory.")
        exit(1)

    # Squeeze out the singleton dimension (Nz, Nx)
    kext = np.squeeze(kext_gaz_raw)  # (Nz, Nx)

    # Convert from dB/km to dB (for 1 km path) – or keep as dB/km for plotting
    # We'll plot dB/km

    # --- Plotting ---
    plt.figure(figsize=(10, 5))
    im = plt.pcolormesh(r / 1000, z / 1000, kext, shading='auto', cmap='viridis', vmin=0)
    plt.colorbar(im, label='Gaseous attenuation (dB/km)')
    plt.xlabel('Range (km)')
    plt.ylabel('Height (km)')
    plt.title(f'W-band ({freq_ghz} GHz) Gaseous Attenuation – Ground Radar Cross-Section')
    plt.ylim(0, 10)
    plt.xlim(0, 50)
    plt.grid(True, color='white', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()