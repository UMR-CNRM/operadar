#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 18:43:26 2026

@author: augros
"""
import numpy as np
import matplotlib.pyplot as plt

def compute_attenuated_zh_3D(kext, zh, model_altitude, radar_altitude):
    # ... (The code you provided in your prompt goes here) ...
    kext = np.asarray(kext)
    zh = np.asarray(zh)
    z = np.asarray(model_altitude)
    Nz, Ny, Nx = zh.shape
    zh_lin = 10 ** (zh / 10.0)
    flip_mask = z[0, ...] > z[-1, ...]
    if np.any(flip_mask):
        z = z.copy(); kext = kext.copy(); zh_lin = zh_lin.copy()
        z[:, flip_mask] = z[::-1, flip_mask]
        kext[:, flip_mask] = kext[::-1, flip_mask]
        zh_lin[:, flip_mask] = zh_lin[::-1, flip_mask]
    dz = np.diff(z, axis=0)
    kext_mean = 0.5 * (kext[1:, ...] + kext[:-1, ...])
    idx = np.abs(z - radar_altitude).argmin(axis=0)
    att = np.zeros_like(zh_lin)
    
    # Note: Your provided loop logic for 'att' has a slight indexing dependency. 
    # For a robust test, we calculate attenuation relative to the radar index.
    for j in range(Ny):
        for i in range(Nx):
            ridx = idx[j,i]
            # Downward from radar
            if ridx > 0:
                segs = kext_mean[:ridx, j, i] * dz[:ridx, j, i]
                att[:ridx, j, i] = np.cumsum(segs[::-1])[::-1]
            # Upward from radar
            if ridx < Nz - 1:
                segs = kext_mean[ridx:, j, i] * dz[ridx:, j, i]
                att[ridx+1:, j, i] = np.cumsum(segs)
                
    zh_att_lin = zh_lin * np.exp(-2.0 * att)
    zh_att = np.full_like(zh, np.nan)
    mask = zh_att_lin > 0
    zh_att[mask] = 10.0 * np.log10(zh_att_lin[mask])
    if np.any(flip_mask):
        zh_att[:, flip_mask] = zh_att[::-1, flip_mask]
    return zh_att

#2D Cross-Section Visualization


import matplotlib.pyplot as plt

# --- (The compute_attenuated_zh_3D function from your prompt is assumed to be defined here) ---

# 1. Setup a 2D Slice (Nx=100, Ny=1 for a cross-section, Nz=100)
Nz, Ny, Nx = 100, 1, 100
x_vec = np.linspace(0, 50, Nx)  # 50 km horizontal distance
z_vec = np.linspace(0, 10000, Nz) # 10 km vertical

# Create 3D arrays (Nz, 1, Nx)
# We add a slight "wave" to the altitude to test the per-column logic
z_3d = np.tile(z_vec[:, None, None], (1, Ny, Nx)) 
for i in range(Nx):
    z_3d[:, 0, i] += 100 * np.sin(i / 10.0) # Terrain/Grid waviness

# Define a cloud layer between 2km and 8km
zh_3d = np.zeros((Nz, Ny, Nx))
kext_3d = np.zeros((Nz, Ny, Nx))

cloud_mask = (z_3d > 2000) & (z_3d < 8000)
zh_3d[cloud_mask] = 35.0      # 35 dBZ cloud
kext_3d[cloud_mask] = 4e-4    # Constant extinction in cloud

# 2. Compute Attenuation for both cases
radar_alt_ground = 0.0
radar_alt_air = 5000.0

zh_att_ground = compute_attenuated_zh_3D(kext_3d, zh_3d, z_3d, radar_alt_ground)
zh_att_air = compute_attenuated_zh_3D(kext_3d, zh_3d, z_3d, radar_alt_air)

# 3. Plotting
fig, axes = plt.subplots(3, 1, figsize=(12, 12), sharex=True)

# Helper to plot the 2D slice
def plot_slice(ax, data, title, radar_alt=None):
    # We take the [:, 0, :] slice to get (Nz, Nx)
    # Using pcolormesh with the specific z_3d coordinates per column
    X, Z = np.meshgrid(x_vec, np.zeros(Nz)) # dummy X for shape
    # Actual coordinates for pcolormesh:
    curr_x = np.tile(x_vec[None, :], (Nz, 1))
    curr_z = z_3d[:, 0, :]
    
    mesh = ax.pcolormesh(curr_x, curr_z, data[:, 0, :], 
                         shading='auto', cmap='pyart_HomeyerRainbow' if 'pyart' in plt.colormaps else 'turbo',
                         vmin=0, vmax=35)
    
    if radar_alt is not None:
        ax.axhline(radar_alt, color='white', linestyle='--', alpha=0.8, label=f'Radar Alt ({radar_alt}m)')
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_ylabel("Altitude (m)")
    return mesh

# Plot 1: Original (Unattenuated)
plot_slice(axes[0], zh_3d, "Original Unattenuated Reflectivity (35 dBZ Cloud)")

# Plot 2: Ground Radar
mesh2 = plot_slice(axes[1], zh_att_ground, "Ground Radar View (Radar at 0m)", radar_alt_ground)

# Plot 3: Aircraft Radar
mesh3 = plot_slice(axes[2], zh_att_air, "Aircraft Radar View (Radar at 5000m)", radar_alt_air)

axes[2].set_xlabel("Horizontal Distance (km)")

# Add a colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
fig.colorbar(mesh3, cax=cbar_ax, label='Reflectivity (dBZ)')

plt.subplots_adjust(right=0.9, hspace=0.3)
plt.show()


# # ==========================================
# # SETUP TEST DATA (3D)
# # ==========================================
# Nz, Ny, Nx = 100, 5, 5
# z_vec = np.linspace(0, 10000, Nz)
# # Create 3D altitude with slight variations per column
# z_3d = np.tile(z_vec[:, None, None], (1, Ny, Nx)) + np.random.uniform(-10, 10, (Nz, Ny, Nx))
# z_3d = np.sort(z_3d, axis=0)

# # Define a cloud from 1km to 9km
# zh_3d = np.full((Nz, Ny, Nx), 30.0) # 30 dBZ
# kext_3d = np.full((Nz, Ny, Nx), 4e-4) # 0.0004 Np/m extinction

# # ==========================================
# # RUN SCENARIOS
# # ==========================================
# # Case 1: Ground Radar (0m)
# zh_ground = compute_attenuated_zh_3D(kext_3d, zh_3d, z_3d, radar_altitude=0.0)

# # Case 2: Aircraft Radar (5000m)
# zh_aircraft = compute_attenuated_zh_3D(kext_3d, zh_3d, z_3d, radar_altitude=5000.0)

# # ==========================================
# # PLOTTING
# # ==========================================
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

# # Select middle column for 1D profile illustration
# yi, xi = Ny//2, Nx//2

# # Plot Ground Radar Scenario
# ax1.plot(zh_3d[:, yi, xi], z_3d[:, yi, xi], 'k--', label='True Reflectivity')
# ax1.plot(zh_ground[:, yi, xi], z_3d[:, yi, xi], 'r-', linewidth=2, label='Attenuated (Ground)')
# ax1.axhline(0, color='blue', linestyle=':', label='Radar Level')
# ax1.set_title("Ground Radar Scenario\n(Radar at 0m)")
# ax1.set_xlabel("Reflectivity (dBZ)")
# ax1.set_ylabel("Altitude (m)")
# ax1.grid(True, alpha=0.3)
# ax1.legend()

# # Plot Aircraft Radar Scenario
# ax2.plot(zh_3d[:, yi, xi], z_3d[:, yi, xi], 'k--', label='True Reflectivity')
# ax2.plot(zh_aircraft[:, yi, xi], z_3d[:, yi, xi], 'g-', linewidth=2, label='Attenuated (Aircraft)')
# ax2.axhline(5000, color='blue', linestyle=':', label='Radar Level (5km)')
# ax2.set_title("Aircraft Radar Scenario\n(Radar at 5000m)")
# ax2.set_xlabel("Reflectivity (dBZ)")
# ax2.grid(True, alpha=0.3)
# ax2.legend()

# plt.tight_layout()
# plt.show()