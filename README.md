# TRAIN

**TRAIN (Toolbox for Reducing Atmospheric InSAR Noise)** is a toolbox for tropospheric correction of InSAR data. The original project was developed by **David P. S. Bekaert and collaborators at the University of Leeds** and provides a collection of MATLAB, shell, and command-line tools for applying atmospheric corrections within InSAR processing workflows, including integration with **StaMPS**.

This repository is a maintained fork of the original TRAIN project:

- Original project: https://github.com/dbekaert/TRAIN
- This fork: https://github.com/limingjia92/TRAIN

The original TRAIN algorithms, workflow design, citation requirements, and GNU GPL license are retained unless explicitly modified for bug fixes or compatibility updates.

## Changes in this fork

The original TRAIN repository has received little development since 2019. This fork modernizes several parts of the codebase for current StaMPS/StaMPS-HPC and ERA5 workflows while preserving the original TRAIN processing convention.

### StaMPS / StaMPS-HPC compatibility

- Improved parameter inheritance from StaMPS `parms.mat`, including wavelength, satellite heading, and acquisition UTC time where available.
- Improved fallback handling for legacy TRAIN parameter files.
- Fixed UTC-time handling used to determine the ERA5 epochs before and after SAR acquisition time.
- Maintained compatibility with the original TRAIN/StaMPS data convention while targeting use with [StaMPS-HPC](https://github.com/limingjia92/StaMPS-HPC).

### ERA5 and current CDS API support

- Added a dedicated `era5_datapath`.
- Updated ERA5 downloading to the current Copernicus Climate Data Store (`cdsapi`) request format.
- Updated ERA5 pressure-level requests for the current CDS NetCDF output.
- Added support for current NetCDF coordinate/dimension names such as `pressure_level` and `valid_time`, while retaining compatibility aliases for older layouts where practical.
- Improved NetCDF dimension handling and pressure-level ordering.
- Added exact ERA5 request manifests (`YYYYMMDDHHMM`) to correctly handle acquisition times crossing UTC midnight.
- Added optimized batch ERA5 downloading to reduce the number of CDS requests.
- Added resumable batch downloads, bounded retries, temporary `.part` files, validation, and atomic output writes.
- Batch files are automatically split back into the TRAIN-compatible layout:

```text
ERA5/
└── YYYYMMDD/
    └── ggapYYYYMMDDHHMM.nc
```

- Existing valid ERA5 files are reused instead of being overwritten automatically.

The updated ERA5 reader and single-epoch CDS download path have been regression-tested using current CDS ERA5 pressure-level NetCDF files.

### GMT / DEM compatibility and bug fixes

- Improved GMT version detection for modern GMT releases.
- Updated `xyz2grd` handling for GMT 5/6 node-count syntax and modern options.
- Corrected DEM binary data-type handling.
- Fixed several historical parameter, file-path, and compatibility issues identified during modernization.

## Installation

Clone this fork directly from GitHub:

```bash
git clone https://github.com/limingjia92/TRAIN.git
cd TRAIN
```

To update an existing clone:

```bash
cd TRAIN
git pull origin master
```

After cloning, configure TRAIN paths and parameters for your local MATLAB/StaMPS environment in the same way as a standard TRAIN installation.

For ERA5 downloading, the Python environment should provide:

```text
cdsapi
xarray
pandas
```

A valid Copernicus Climate Data Store account and `~/.cdsapirc` configuration are also required.

## Citation

If you use TRAIN, please cite the original TRAIN publication:

> Bekaert, D. P. S., Walters, R. J., Wright, T. J., Hooper, A. J., and Parker, D. J. (2015), Statistical comparison of InSAR tropospheric correction techniques, *Remote Sensing of Environment*. https://doi.org/10.1016/j.rse.2015.08.035

When using specific correction methods, please also cite the corresponding original methodological publications. For example, for the spatially variable power-law method:

> Bekaert, D. P. S., Hooper, A. J., and Wright, T. J. (2015), A spatially-variable power-law tropospheric correction technique for InSAR data, *Journal of Geophysical Research: Solid Earth*. https://doi.org/10.1029/2014JB011558

Please also acknowledge the original TRAIN project when using this fork.

## License

TRAIN is distributed under the **GNU General Public License (GPL)**. See the repository license files for details.

## Acknowledgements

The original TRAIN development and community contributions are credited to David P. S. Bekaert and collaborators, including Richard J. Walters, Hannes Bathke, Simran Sangha, Tim J. Wright, Andy J. Hooper, Doug J. Parker, Zhenhong Li, the Leeds InSAR group, COMET members, and the wider TRAIN community.

This fork builds on that work and focuses on maintaining compatibility with modern StaMPS/StaMPS-HPC, ERA5/CDS, MATLAB, NetCDF, and GMT workflows.
::: ​​