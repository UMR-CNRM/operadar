"""Access to the lookup tables for the plot tools (NetCDF <table>.nc, or text tables not converted).

Uses operadar/read/table_io.py, imported from the installed operadar package or from this repository."""
import os
import sys

try:
    from operadar.read.table_io import table_path, read_table_params, read_table_dataframe, open_table
except ImportError:  # operadar not installed: use the module of this repository
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "operadar", "read"))
    from table_io import table_path, read_table_params, read_table_dataframe, open_table


def table_exists(path: str) -> bool:
    """True if the table exists as <path>.nc or as the text file <path>."""
    return os.path.isfile(path) or os.path.isfile(path + ".nc")
