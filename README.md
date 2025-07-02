# GRAF Reforecast Documentation

## Introduction

This repository provides documentation and tools for accessing and visualizing GRAF (Global Reforecasts from MPAS - Model for Prediction Across Scales) reforecast data. GRAF is a weather forecasting system that features mesh refinement over the US and Europe, allowing for higher resolution forecasts in these regions.

The reforecast dataset contains historical forecasts generated using the GRAF model, which can be valuable for research, verification, and calibration purposes.

## Repository Contents

- `demo_read_zarr_s3.py`: Main demo script for reading and visualizing GRAF reforecast data from S3 storage
- `GRAF_reforecast_initial_condition_dates_v2.txt`: List of available initial condition dates for the reforecasts
- `dateutils.py`: Utility functions for date manipulation
- `rpm4km.static.latlon.nc`: NetCDF file containing lat/lon coordinates for the unstructured mesh grid points
- `requirements.txt`: List of Python dependencies required to run the demo

The GRAF reforecast data is stored in Zarr format in an S3 data store, which allows for efficient access to specific subsets of the data.

## Prerequisites

- Python 3.12.0 or any compatible python version

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/your-username/graf-reforecast-documentation.git
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

The demo script `demo_read_zarr_s3.py` allows you to read and visualize specific GRAF reforecast data. The script requires three command-line arguments:

```
python demo_read_zarr_s3.py <initial_condition_date> <lead_time> <variable>
```

Where:
- `<initial_condition_date>`: Initial condition date in YYYYMMDDHH format
- `<lead_time>`: Lead time in hours (can include fractions for sub-hourly data)
- `<variable>`: Variable name to plot

### Example

```
python demo_read_zarr_s3.py 2023122012 12 apcp_bucket
```

This command will:
1. Access the GRAF reforecast data for the initial condition date December 20, 2023 at 12Z
2. Extract the total precipitation field at the 12-hour forecast lead time
3. Create a plot of the data and save it as a PNG file

The script will output the path to the generated plot file upon completion.
