#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4D Hypercube Interpolation for Scattering Coefficients.

This module provides multilinear interpolation in 4D space (Elevation, Temperature, 
P3-parameter, Content) for radar scattering simulations.

Created on Thu Apr 23 11:13:59 2026
@author: augrosc & davidcl
"""

import numpy as np
from typing import List, Tuple, Dict, Union
import warnings


def get_bounds(val: np.ndarray, val_min: float, val_max: float, step: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute lower and upper grid bounds for interpolation.
    
    For values exactly on the maximum boundary, ensures the interval
    [val_max - step, val_max] is used to avoid division by zero.
    
    Parameters
    ----------
    val : np.ndarray
        Input values (can be any shape)
    val_min, val_max : float
        Grid boundaries
    step : float
        Grid step size
        
    Returns
    -------
    val_inf, val_sup : Tuple[np.ndarray, np.ndarray]
        Lower and upper bound arrays (same shape as val)
    """
    # Clip values to grid range first
    val_clipped = np.clip(val, val_min, val_max)
    
    # Calculate normalized coordinate
    norm = (val_clipped - val_min) / step
    idx = np.floor(norm).astype(int)
    
    # For values at the maximum boundary, use the last valid interval
    n_intervals = int(np.round((val_max - val_min) / step))
    idx = np.clip(idx, 0, n_intervals - 1)
    
    val_inf = idx * step + val_min
    val_sup = val_inf + step
    
    # Ensure sup doesn't exceed max (numerical safety)
    val_sup = np.minimum(val_sup, val_max)
    
    return val_inf, val_sup


def get_linear_index(T: np.ndarray, E: np.ndarray, P: np.ndarray, M: np.ndarray,
                     T_min: float, E_min: float, P_min: float, M_min: float,
                     step_T: float, step_E: float, step_P: float, step_M: float,
                     n_T: int, n_E: int, n_P: int, n_M: int) -> np.ndarray:
    """
    Compute linear index into a flattened 4D grid (T, E, P, M) with C-order.
    
    Indices are clipped to valid range [0, n-1] for safety.
    """
    # Calculate integer indices
    iT = np.floor((T - T_min) / step_T).astype(int)
    iE = np.floor((E - E_min) / step_E).astype(int)
    iP = np.floor((P - P_min) / step_P).astype(int)
    iM = np.floor((M - M_min) / step_M).astype(int)

    # Clip to valid grid indices
    iT = np.clip(iT, 0, n_T - 1)
    iE = np.clip(iE, 0, n_E - 1)
    iP = np.clip(iP, 0, n_P - 1)
    iM = np.clip(iM, 0, n_M - 1)

    # Convert to linear index: T is slowest varying, M is fastest
    return ((iT * n_E + iE) * n_P + iP) * n_M + iM


def hypercube_interpolation(
    tableDict: Dict[str, Union[Dict, np.ndarray]],
    hydrometeor: str,
    colName: str,
    elev: np.ndarray,
    T: np.ndarray,
    P3: np.ndarray,
    content: np.ndarray,
    columns_to_retrieve: List[str],
    test_mode: bool = False,
) -> Dict[str, np.ndarray]:
    """
    Perform 4D multilinear interpolation on scattering coefficient tables.
    
    Interpolates in the space: (Elevation, Temperature, P3, Content) where
    Content is always log10-transformed, and P3 is either linear (Fw) or 
    log-transformed (Nc).
    
    Parameters
    ----------
    tableDict : dict
        Dictionary containing lookup tables with keys like:
        'Tcmin', 'Tcmax', 'Tcstep', 'ELEVmin', etc. for grid geometry
        and scattering coefficient arrays per hydrometeor.
    hydrometeor : str
        Hydrometeor type (key for tableDict sub-dictionaries)
    colName : str
        Column type for P3: 'Fw' (linear) or 'Nc' (log10)
    elev : np.ndarray
        Elevation angles [degrees]
    T : np.ndarray
        Temperature [°C]
    P3 : np.ndarray
        Third parameter (Fw or Nc depending on colName)
    content : np.ndarray
        Water content [kg/m^3 or similar]
    columns_to_retrieve : List[str]
        List of scattering coefficient names to interpolate (e.g., ['sighh', 'sighv'])
    test_mode : bool, optional
        If True, prints debug information for random sample points
        
    Returns
    -------
    dict
        Dictionary with interpolated values for each requested column
        
    Raises
    ------
    ValueError
        If input arrays have mismatched shapes or invalid colName
    KeyError
        If required keys are missing from tableDict
    """
    # --- 1. Input Validation ---
    if colName not in ['Fw', 'Nc']:
        raise ValueError(f"colName must be 'Fw' or 'Nc', got {colName}")
    
    # Ensure all input arrays have the same shape
    input_arrays = {'elev': elev, 'T': T, 'P3': P3, 'content': content}
    shapes = [arr.shape for arr in input_arrays.values()]
    if len(set(shapes)) != 1:
        raise ValueError(f"Input arrays must have same shape. Got shapes: {dict(zip(input_arrays.keys(), shapes))}")
    
    # Check required tableDict keys exist
    required_keys = ['Tcmin', 'Tcmax', 'Tcstep', 'ELEVmin', 'ELEVmax', 'ELEVstep',
                     'expMmin', 'expMmax', 'expMstep']
    for key in required_keys:
        if key not in tableDict or hydrometeor not in tableDict[key]:
            raise KeyError(f"Missing required key: tableDict['{key}']['{hydrometeor}']")
    
    # --- 2. Variable Preparation ---
    # Log-transform content (mass)
    expntM = np.full_like(content, -100.0)  # Default for zero content
    mask_positive = content > 0
    expntM[mask_positive] = np.log10(content[mask_positive])
    
    # Handle P3 variable based on type
    if colName == 'Fw':
        P_min_grid = tableDict['Fwmin'][hydrometeor]
        P_max_grid = tableDict['Fwmax'][hydrometeor]
        step_P = tableDict['Fwstep'][hydrometeor]
        P = np.copy(P3)  # Keep linear
    else:  # 'Nc'
        P_min_grid = tableDict['expCCmin'][hydrometeor]
        P_max_grid = tableDict['expCCmax'][hydrometeor]
        step_P = tableDict['expCCstep'][hydrometeor]
        # Log-transform, handling zeros/negatives
        P = np.full_like(P3, P_min_grid)
        mask_positive = P3 > 0
        P[mask_positive] = np.log10(P3[mask_positive])
    
    # --- 3. Grid Geometry ---
    T_min, T_max, step_T = (tableDict['Tcmin'][hydrometeor], 
                            tableDict['Tcmax'][hydrometeor], 
                            tableDict['Tcstep'][hydrometeor])
    elev_min, elev_max, step_elev = (tableDict['ELEVmin'][hydrometeor], 
                                     tableDict['ELEVmax'][hydrometeor], 
                                     tableDict['ELEVstep'][hydrometeor])
    expntM_min, expntM_max, step_expntM = (tableDict['expMmin'][hydrometeor], 
                                           tableDict['expMmax'][hydrometeor], 
                                           tableDict['expMstep'][hydrometeor])
    
    # Calculate grid dimensions
    n_T = int(np.round((T_max - T_min) / step_T)) + 1
    n_elev = int(np.round((elev_max - elev_min) / step_elev)) + 1
    n_P = int(np.round((P_max_grid - P_min_grid) / step_P)) + 1
    n_expntM = int(np.round((expntM_max - expntM_min) / step_expntM)) + 1
    
    # --- 4. Compute Interpolation Bounds ---
    T_inf, T_sup = get_bounds(T, T_min, T_max, step_T)
    elev_inf, elev_sup = get_bounds(elev, elev_min, elev_max, step_elev)
    P_inf, P_sup = get_bounds(P, P_min_grid, P_max_grid, step_P)
    expntM_inf, expntM_sup = get_bounds(expntM, expntM_min, expntM_max, step_expntM)
    
    # --- 5. Compute Interpolation Weights ---
    # Avoid division by zero (can happen if step is 0 or values identical)
    eps = 1e-12
    
    wT = np.clip((T - T_inf) / np.maximum(T_sup - T_inf, eps), 0, 1)
    wE = np.clip((elev - elev_inf) / np.maximum(elev_sup - elev_inf, eps), 0, 1)
    wP = np.clip((P - P_inf) / np.maximum(P_sup - P_inf, eps), 0, 1)
    wM = np.clip((expntM - expntM_inf) / np.maximum(expntM_sup - expntM_inf, eps), 0, 1)
    
    # --- 6. Generate 16 Corners of 4D Hypercube ---
    # Each corner is a combination of (lower/upper) for each dimension
    # 0 = lower (inf), 1 = upper (sup)
    corners = np.array([[i, j, k, l] 
                       for i in [0, 1] 
                       for j in [0, 1] 
                       for k in [0, 1] 
                       for l in [0, 1]])  # Shape: (16, 4)
    
    # --- 7. Compute Corner Coordinates and Weights ---
    # Vectorized computation for all 16 corners
    # Shape: (16, N) where N is number of points
    # N = elev.size
    
    # # Expand dims for broadcasting: (16, 1) vs (N,)
    # corner_T = np.where(corners[:, 0:1] == 0, T_inf.flat, T_sup.flat)  # (16, N)
    # corner_E = np.where(corners[:, 1:2] == 0, elev_inf.flat, elev_sup.flat)
    # corner_P = np.where(corners[:, 2:3] == 0, P_inf.flat, P_sup.flat)
    # corner_M = np.where(corners[:, 3:4] == 0, expntM_inf.flat, expntM_sup.flat)
    
    N = elev.size
    T_inf_r = T_inf.ravel()
    T_sup_r = T_sup.ravel()
    elev_inf_r = elev_inf.ravel()
    elev_sup_r = elev_sup.ravel()
    P_inf_r = P_inf.ravel()
    P_sup_r = P_sup.ravel()
    expntM_inf_r = expntM_inf.ravel()
    expntM_sup_r = expntM_sup.ravel()
    
    # Expand dims for broadcasting: (16, 1) vs (N,)
    corner_T = np.where(corners[:, 0:1] == 0, T_inf_r, T_sup_r)
    corner_E = np.where(corners[:, 1:2] == 0, elev_inf_r, elev_sup_r)
    corner_P = np.where(corners[:, 2:3] == 0, P_inf_r, P_sup_r)
    corner_M = np.where(corners[:, 3:4] == 0, expntM_inf_r, expntM_sup_r)
    
    # Weights for each corner: product of individual dimension weights
    # (1-w) for lower, w for upper
    weight_T = np.where(corners[:, 0:1] == 0, 1.0 - wT[None, :], wT[None, :])
    weight_E = np.where(corners[:, 1:2] == 0, 1.0 - wE[None, :], wE[None, :])
    weight_P = np.where(corners[:, 2:3] == 0, 1.0 - wP[None, :], wP[None, :])
    weight_M = np.where(corners[:, 3:4] == 0, 1.0 - wM[None, :], wM[None, :])
    
    corner_weights = weight_T * weight_E * weight_P * weight_M  # (16, N)
    
    # --- 8. Compute Linear Indices ---
    idx_all = get_linear_index(
        corner_T, corner_E, corner_P, corner_M,
        T_min, elev_min, P_min_grid, expntM_min,
        step_T, step_elev, step_P, step_expntM,
        n_T, n_elev, n_P, n_expntM
    )
    
    # --- 9. Interpolate Each Requested Column ---
    scatCoefsDict = {}
    
    for column in columns_to_retrieve:
        if column not in tableDict:
            warnings.warn(f"Column '{column}' not found in tableDict, returning zeros")
            scatCoefsDict[column] = np.zeros_like(elev)
            continue
            
        if hydrometeor not in tableDict[column]:
            warnings.warn(f"Hydrometeor '{hydrometeor}' not found in tableDict['{column}'], returning zeros")
            scatCoefsDict[column] = np.zeros_like(elev)
            continue
        
        table_values = tableDict[column][hydrometeor]
        
        # Fetch values at all 16 corners: shape (16, N)
        vals_at_corners = table_values[idx_all]
        
        # Weighted sum: sum over 16 corners
        interpolated = np.sum(corner_weights * vals_at_corners, axis=0)
        
        # Reshape back to original input shape
        scatCoefsDict[column] = interpolated.reshape(elev.shape)
    
    # --- 10. Debug Output ---
    if test_mode:
        np.random.seed(42)  # Reproducible sampling
        sample_size = min(5, N)
        sample_indices = np.random.choice(N, size=sample_size, replace=False)
        
        print(f"\n=== Debug Info (showing {sample_size} random points) ===")
        print(f"Grid dims: T={n_T}, E={n_elev}, P={n_P}, M={n_expntM}")
        
        for i, idx in enumerate(sample_indices):
            print(f"\nPoint {i+1} (flat idx {idx}):")
            print(f"  Inputs: elev={elev.flat[idx]:.2f}, T={T.flat[idx]:.2f}, "
                  f"P={P.flat[idx]:.4f}, M={expntM.flat[idx]:.4f}")
            for col in columns_to_retrieve[:3]:  # Limit output
                if col in scatCoefsDict:
                    print(f"  {col} = {scatCoefsDict[col].flat[idx]:.6e}")
    
    return scatCoefsDict


def create_test_table() -> Dict:
    """
    Create a synthetic 4D lookup table for testing.

    For Fw mode, interpolation axis uses linear P3 in [0,1].
    For Nc mode, interpolation axis uses log10(P3) in [expCCmin, expCCmax].

    We therefore store two different coefficients:
      - 'test_val_fw' valid for colName='Fw'
      - 'test_val_nc' valid for colName='Nc'
    """
    # Grid parameters (same for T, E, M)
    T_min, T_max, step_T = -20.0, 20.0, 10.0
    E_min, E_max, step_E = 0.0, 90.0, 30.0

    # Linear P axis for Fw
    P_fw_min, P_fw_max, step_P_fw = 0.0, 1.0, 0.25

    # Log10(P) axis for Nc
    # Must match what tests set for expCCmin/expCCmax/expCCstep
    expCC_min, expCC_max, step_P_nc = 0.0, 3.0, 1.0  # log10(P3)

    # M axis is log10(content), so grid uses "log content"
    M_min, M_max, step_M = -3.0, 1.0, 1.0

    n_T = int((T_max - T_min) / step_T) + 1
    n_E = int((E_max - E_min) / step_E) + 1
    n_P_fw = int((P_fw_max - P_fw_min) / step_P_fw) + 1
    n_P_nc = int((expCC_max - expCC_min) / step_P_nc) + 1
    n_M = int((M_max - M_min) / step_M) + 1

    T_grid = np.linspace(T_min, T_max, n_T)
    E_grid = np.linspace(E_min, E_max, n_E)
    P_fw_grid = np.linspace(P_fw_min, P_fw_max, n_P_fw)
    P_nc_grid = np.linspace(expCC_min, expCC_max, n_P_nc)
    M_grid = np.linspace(M_min, M_max, n_M)

    # Build Fw mesh
    TT, EE, PP_fw, MM = np.meshgrid(T_grid, E_grid, P_fw_grid, M_grid, indexing='ij')
    test_data_fw = TT + 2.0 * EE + 3.0 * PP_fw + 4.0 * MM
    test_data_fw_flat = test_data_fw.flatten()

    # Build Nc mesh (PP is log10(P3))
    TT, EE, PP_nc, MM = np.meshgrid(T_grid, E_grid, P_nc_grid, M_grid, indexing='ij')
    test_data_nc = TT + 2.0 * EE + 3.0 * PP_nc + 4.0 * MM
    test_data_nc_flat = test_data_nc.flatten()

    return {
        # Grid geometry
        'Tcmin': {'rain': T_min}, 'Tcmax': {'rain': T_max}, 'Tcstep': {'rain': step_T},
        'ELEVmin': {'rain': E_min}, 'ELEVmax': {'rain': E_max}, 'ELEVstep': {'rain': step_E},
        'expMmin': {'rain': M_min}, 'expMmax': {'rain': M_max}, 'expMstep': {'rain': step_M},

        # P axes
        'Fwmin': {'rain': P_fw_min}, 'Fwmax': {'rain': P_fw_max}, 'Fwstep': {'rain': step_P_fw},
        'expCCmin': {'rain': expCC_min}, 'expCCmax': {'rain': expCC_max}, 'expCCstep': {'rain': step_P_nc},

        # Coefficients (mode-specific)
        'test_val_fw': {'rain': test_data_fw_flat.copy()},
        'test_val_nc': {'rain': test_data_nc_flat.copy()},
    }


def run_tests():
    """Comprehensive test suite for hypercube_interpolation."""
    print("Running hypercube_interpolation tests...")
    
    # Create test data
    table = create_test_table()
    
    print("Columns in synthetic table:", 
          [k for k in table.keys() if isinstance(table[k], dict)])
    
    # Test 1: Exact interpolation at grid points
    print("\n1. Testing exact interpolation at grid points...")
    # Pick specific grid point: T=0, E=30, P=0.5, M=-1
    elev_test = np.array([30.0])
    T_test = np.array([0.0])
    P3_test = np.array([0.5])
    content_test = np.array([10**(-1.0)])  # M = -1
    
    result = hypercube_interpolation(
        table, 'rain', 'Fw', elev_test, T_test, P3_test, content_test,
        columns_to_retrieve=['test_val_fw']
    )
    
    expected = 0.0 + 2.0*30.0 + 3.0*0.5 + 4.0*(-1.0)  # 60 + 1.5 - 4 = 57.5
    actual = result['test_val_fw'][0]
    
    assert np.allclose(actual, expected, rtol=1e-10), \
        f"Grid point test failed: expected {expected}, got {actual}"
    print(f"   ✓ Grid point test passed (value={actual:.6f})")
    
    # Test 2: Midpoint interpolation (should be exact for linear function)
    print("\n2. Testing midpoint interpolation...")
    # Midway between grid points: T=5 (between 0 and 10), E=15 (between 0 and 30), etc.
    elev_test = np.array([15.0])
    T_test = np.array([5.0])
    P3_test = np.array([0.375])  # Between 0.25 and 0.5
    content_test = np.array([10**(-2.0)])  # M = -2 (between -3 and -1)
    
    result = hypercube_interpolation(
        table, 'rain', 'Fw', elev_test, T_test, P3_test, content_test,
        columns_to_retrieve=['test_val_fw']
    )
    
    expected = 5.0 + 2.0*15.0 + 3.0*0.375 + 4.0*(-2.0)  # 5 + 30 + 1.125 - 8 = 28.125
    actual = result['test_val_fw'][0]
    
    assert np.allclose(actual, expected, rtol=1e-10), \
        f"Midpoint test failed: expected {expected}, got {actual}"
    print(f"   ✓ Midpoint test passed (value={actual:.6f})")
    
    # Test 3: Boundary handling (extrapolation should clamp)
    print("\n3. Testing boundary clamping...")
    # Values outside grid
    elev_test = np.array([100.0])  # Above max (90)
    T_test = np.array([-30.0])     # Below min (-20)
    P3_test = np.array([0.5])
    content_test = np.array([10**(-1.0)])
    
    result = hypercube_interpolation(
        table, 'rain', 'Fw', elev_test, T_test, P3_test, content_test,
        columns_to_retrieve=['test_val_fw']
    )
    
    # Should use values at T=-20, E=90
    expected = (-20.0) + 2.0*90.0 + 3.0*0.5 + 4.0*(-1.0)  # -20 + 180 + 1.5 - 4 = 157.5
    actual = result['test_val_fw'][0]
    
    assert np.allclose(actual, expected, rtol=1e-10), \
        f"Boundary test failed: expected {expected}, got {actual}"
    print(f"   ✓ Boundary clamping test passed (value={actual:.6f})")
    
    # Test 4: Nc mode (log-transformed P3)
    print("\n4. Testing 'Nc' mode (log P3)...")
    
    elev_test = np.array([0.0])
    T_test = np.array([-20.0])
    P3_test = np.array([100.0])          # log10(100)=2.0
    content_test = np.array([10**(-3.0)])  # M=-3
    
    result = hypercube_interpolation(
        table, 'rain', 'Nc',
        elev_test, T_test, P3_test, content_test,
        columns_to_retrieve=['test_val_nc']
    )
    
    expected = (-20.0) + 2.0*0.0 + 3.0*2.0 + 4.0*(-3.0)
    actual = result['test_val_nc'][0]
    
    assert np.allclose(actual, expected, rtol=1e-10), \
        f"Nc mode test failed: expected {expected}, got {actual}"
    print(f"   ✓ Nc mode test passed (value={actual:.6f})")
    
    # Test 5: Array input consistency (linear function should match exactly)
    print("\n5. Testing array inputs (linear function exactness)...")
    np.random.seed(123)
    n_points = 1000
    elev_arr = np.random.uniform(0, 90, n_points)
    T_arr = np.random.uniform(-20, 20, n_points)
    P3_arr = np.random.uniform(0, 1, n_points)          # Fw uses linear P3
    content_arr = 10**np.random.uniform(-3, 1, n_points)  # M = log10(content)
    
    result = hypercube_interpolation(
        table, 'rain', 'Fw',
        elev_arr, T_arr, P3_arr, content_arr,
        columns_to_retrieve=['test_val_fw']
    )
    
    # Expected: TT + 2*EE + 3*PP + 4*MM
    # where MM = log10(content_arr)
    M_expected = np.log10(content_arr)
    expected = T_arr + 2.0 * elev_arr + 3.0 * P3_arr + 4.0 * M_expected
    
    assert result['test_val_fw'].shape == (n_points,)
    assert np.allclose(result['test_val_fw'], expected, rtol=1e-10, atol=1e-12), \
        "Array interpolation consistency check failed (linear exactness)"
    print(f"   ✓ Array input test passed ({n_points} points)")
    
    # Test 6: Zero content handling
    print("\n6. Testing zero content handling...")
    content_zero = np.array([0.0, 1e-12, 1e-6])  # Zero and near-zero
    elev_test = np.array([0.0, 0.0, 0.0])
    T_test = np.array([0.0, 0.0, 0.0])
    P3_test = np.array([0.5, 0.5, 0.5])
    
    result = hypercube_interpolation(
        table, 'rain', 'Fw', elev_test, T_test, P3_test, content_zero,
        columns_to_retrieve=['test_val']
    )
    # Should not crash, should return values at M_min (-3)
    print(f"   ✓ Zero content handled (values: {result['test_val']})")
    
    # Test 7: Error handling
    print("\n7. Testing error handling...")
    try:
        hypercube_interpolation(
            table, 'rain', 'InvalidMode', elev_test, T_test, P3_test, content_zero,
            columns_to_retrieve=['test_val']
        )
        assert False, "Should have raised ValueError"
    except ValueError as e:
        print(f"   ✓ Invalid colName caught: {e}")
    
    try:
        bad_elev = np.array([1, 2])  # Wrong shape
        hypercube_interpolation(
            table, 'rain', 'Fw', bad_elev, T_test, P3_test, content_zero,
            columns_to_retrieve=['test_val']
        )
        assert False, "Should have raised ValueError"
    except ValueError as e:
        print(f"   ✓ Shape mismatch caught")
    
    print("\n=== All tests passed! ===")


if __name__ == "__main__":
    run_tests()
    
    # Example usage demonstration
    print("\n" + "="*50)
    print("Example Usage:")
    print("="*50)
    
    table = create_test_table()
    
    # Single point example
    result = hypercube_interpolation(
        tableDict=table,
        hydrometeor='rain',
        colName='Fw',
        elev=np.array([45.0]),
        T=np.array([10.0]),
        P3=np.array([0.5]),
        content=np.array([0.1]),  # 0.1 kg/m^3
        columns_to_retrieve=['sighh', 'sighv'],
        test_mode=True
    )
    
    print(f"\nInterpolated sighh: {result['sighh'][0]:.4e}")
    print(f"Interpolated sighv: {result['sighv'][0]:.4e}")