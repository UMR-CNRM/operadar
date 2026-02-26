from setuptools import setup, find_packages

setup(
    name='operadar',
    version='0.0.1',
    author='C. Augros, C. David',
    description='Radar forward operator developed at the National Centre for Meteorological Research, France',
    url='https://github.com/UMR-CNRM/operadar.git',
    packages=find_packages(exclude=("tests", "plot_tools", "docs")),
    python_requires='>=3.11, <3.13',
    install_requires=[
        "numpy",
        "xarray",
        "pandas",
        "netCDF4",
        "matplotlib",
        "dask",
        # cartopy and epygram should be installed via conda or pip respectively.
    ],
    include_package_data=True,
)
