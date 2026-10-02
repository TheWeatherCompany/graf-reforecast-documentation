# GRAF Reforecast Documentation

## Introduction

This repository provides documentation and example code for accessing and
visualizing the GRAF reforecast data set. GRAF (Global high-Resolution
Atmospheric Forecasting) is The Weather Company's version of NCAR's Model
for Prediction Across Scales (MPAS), run globally with mesh refinement to
approximately 4 km over the US, the Caribbean Basin, and Europe (15 km
elsewhere). The reforecast comprises 1,836 cases (roughly 5 years of
output) generated from ERA5-initialized runs for dates spanning 2004-2024,
selected to emphasize high-impact weather.

The data set is described in full in:

> Hamill, T. M., A. Subramaniam, M. Flora, R. R. Prasanna Kumar, K. Kashinath,
> C. Ponder, T. Ge, S. Khajehie, J. Wong, B. Wilt, and P. Neilley, 2026:
> *reGRAF: a global MPAS reforecast data set with convection-allowing
> refinements over the US and Europe.* [arXiv:2609.40133](https://arxiv.org/abs/2609.40133)

The data are hosted at no cost through AWS's Open Data program; see the
[registry page](https://registry.opendata.aws/graf-reforecast/) and
[AWS Marketplace listing](https://aws.amazon.com/marketplace/pp/prodview-j422yaekib6ya)
for access details.

## Recent updates to the data set (October 2026)

The reforecast has been reprocessed into a corrected "canonical" copy,
stored under `s3://twc-graf-reforecast/2026-05-01.canonical/`. This fixes
several issues present in the original release:

- **Temporal ordering.** The original zarr stores did not always have
  their `Time` coordinate monotonically increasing with lead time; reading
  a timestep by index could return the wrong valid time. The canonical
  copy is written in lead-time order, and the example scripts in this
  repository additionally decode and sort the `Time` coordinate
  themselves as a safeguard.
- **Compression.** The original data were stored uncompressed. The
  canonical copy uses an error-bounded lossy compression pipeline: a
  z-score Normalize filter, [BitRound](https://doi.org/10.1038/s43588-021-00156-2)
  truncation to 10 mantissa bits, a byte-shuffle filter, and Zlib/DEFLATE
  level 9. This reduces storage volume by roughly a factor of 6 while
  keeping more precision than is used by the reduced-precision (TF32 /
  float16 / bfloat16) formats typical of downstream deep-learning
  training pipelines. For precipitation and other variables where exact
  zero values and domain integrals matter, the Normalize filter's mean is
  fixed at zero rather than the variable's z-score mean, which avoids
  bit-truncation zero-drift.
- **Missing cases.** The original release was missing a number of case
  dates. The canonical copy fills these in, bringing the archive to its
  full complement of 1,836 cases.

See the arXiv preprint above (Sections 3-4) for full details on case
selection and data storage, and Section 5 for a known issue with
spurious small-scale vertical velocities that develop in some reforecasts
over long lead times.

## Repository Contents

- `demo_read_zarr_s3.py`: Demo script for reading and plotting a single
  variable/lead-time from the canonical GRAF reforecast data on S3
- `demo_animate_zarr_s3.py`: Demo script that loops over every available
  lead time for one case/variable and assembles the frames into an
  animated GIF
- `GRAF_reforecast_initial_condition_dates_v2.txt`: List of available
  initial condition dates for the reforecast, and each case's maximum
  lead time
- `dateutils.py`: Utility functions for date manipulation
- `rpm4km.static.latlon.nc`: NetCDF file containing lat/lon coordinates
  for the unstructured MPAS mesh grid points
- `requirements.txt`: Python dependencies required to run the demos

The GRAF reforecast data are stored in zarr format in an S3 data store,
which allows for efficient access to specific subsets of the data.

## Prerequisites

- Python 3.12.0 or any compatible python version

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/TheWeatherCompany/graf-reforecast-documentation.git
   cd graf-reforecast-documentation
   ```

2. Set up a Python environment with the required version:
   ```
   # If using pyenv
   pyenv install 3.12.0
   pyenv local 3.12.0
   ```

3. Install the required dependencies:
   ```
   pip install -r requirements.txt
   ```

## Usage

### Plotting a single variable/lead time

```
python demo_read_zarr_s3.py <initial_condition_date> <lead_time> <variable>
```

Where:
- `<initial_condition_date>`: Initial condition date in YYYYMMDDHH format
- `<lead_time>`: Lead time in hours (can include fractions for sub-hourly data)
- `<variable>`: Variable name to plot (run the script with an invalid name
  to print the full list)

Example:

```
python demo_read_zarr_s3.py 2023122012 12 apcp_bucket
```

This command will:
1. Access the GRAF reforecast data for the initial condition date December 20, 2023 at 12Z
2. Extract the total precipitation field at the 12-hour forecast lead time
3. Create a plot of the data and save it as a PNG file

The script will output the path to the generated plot file upon completion.

### Animating a whole forecast

```
python demo_animate_zarr_s3.py <initial_condition_date> <variable>
```

Example:

```
python demo_animate_zarr_s3.py 2004010112 t2m
```

This reads every available lead time for that case/variable, plots each
one, and assembles the frames into an animated GIF (`GRAFR_<variable>_<initial_condition_date>_animation.gif`).
