#This script is to compare two Tmatrix_tables.
import numpy as np
import pandas as pd


exp_name1 = "default"
exp_name2 = "default_copied"
hydros = ["rr", "ss", "gg", "cl", "ii", "wg"]
band = "K"

exp_physics1 = 'ICE3'
exp_physics2 = 'ICE3'

# OPEN TmatCoefDiff
tables1 = {}
tables2 = {}

print(f"reading table {exp_name1} and {exp_name2}")
for h in hydros:
    tables1[h] = pd.read_csv(
        f"../../tables_generator/tables/{exp_name1}/TmatCoefDiff_{band}{h}",
        sep=r"\s+",
        skiprows=2,)
    
    tables2[h] = pd.read_csv(
        f"../../tables_generator/tables/{exp_name2}/TmatCoefDiff_{band}{h}",
        sep=r"\s+",
        skiprows=2,)

print(f"Calculating differences for TmatCoefDiff_{band}")
diffs = {}
for h in hydros:
    diffs[h] = (tables1[h] - tables2[h]).abs()

    max_diff = np.nanmax(diffs[h].to_numpy())
    print(f"{h.upper()}: max diff = {max_diff}")