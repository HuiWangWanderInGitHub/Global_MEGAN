# Global Version Model of Emissions of Gases and Aerosols from Nature (Global-MEGAN)

> `megan` is a model used to estimate biogenic VOCs emissions, it has been developing by Alex Guenther and his group. This is a new version of MEGAN developed by Dr. Hui Wang for global biogenic reactive gases simulations based on Ram's code.

## Dependencies:

+  Fortran compiler
+  NetCDF library
+  MPI library

> &#9888; At this moment it only works on UNIX/Linux O.S. 

## Get MEGAN input data:

In order to run megan you will need meteorological and land data. 

Data required:
+ ERA5 NetCDF file. <!-- with the following variables: 'XLAT', 'XLONG', 'Times', 'MAPFAC_M', 'ISLTYP', 'U10', 'V10', 'T2', 'SWDOWN', 'PSFC', 'Q2', 'RAINNC', 'LAI', 'SMOIS', 'TSLB'. -->
+ Land vegetation fields needed. This version is driven by the gridded emission factors based on the vegetation species and corresponding emission factors. 

## Build
Go to the `src` directory:

`> cd src/`

Edit the Makefile to set the compiler and path to NetCDF lib and include files.

`> make`

If the compilation is successful, the executable `meganv3.3.exe` should be created at the `./exe` directory.

> &#9888; Note that the variables must be adjusted to match the appropriate values for your system. Check your `nc-config --libdir` and `nc-config --includedir` to fill NetCDF flags on Makefile.

## Run

Edit the namelist `namelist_megan` that contains the following variables:


```fortran
&megan_nl
  start_date='__START__',   ! YYYY-MM-DD HH:MM:SS
  end_date  ='__END__',   ! YYYY-MM-DD HH:MM:SS
  output_frequency = '24h', ! Mean over 1h, 3h, 6h, 24h or monthly; one file per month

  met_file_path = '/glade/campaign/univ/ucir0060/ERA5_land/unzip_era5land_2000-2020'
  lai_file_path = '/glade/campaign/univ/ucir0060/LAI4g_Regridded_0.1deg'
  pft_file_path = '/glade/campaign/univ/ucir0060/regridded_GF' 
  ef_file_path = '/glade/campaign/univ/ucir0060/Regridded_0.1deg_EF' 
  output_path   = '/glade/campaign/univ/ucir0060/Global_megan_output_2020' 
  co2_file = "co2_mm_gl.csv"

  met_files ='era5land_0.1degree_<time>.nc',            ! Meteo file path
  pft_files ='ESACCI_GF_0.1_Degree_<time>.nc'           ! pft file path
  ef_file  ='Gridded_GR_EF.nc'                          ! gridded emission factor file
  lai_files ='LAI4g_v1.2_0.1_Degree_<time>.nc'          ! LAI file path

  nlai= '24',                                           ! LAI record number
  lai_scale_factor = 0.1,                               ! LAI file scale factor

  run_flower = .false.                                  ! run flower emission (5% of the leaf level)
  run_litter = .false.                                  ! run litter emission (5% of the leaf level)
  real_spinup_met = .true.                              ! spin up for the long-term meterology
  output_ef_file  = .false.                             ! write the gridded emission factor
  run_co2 = .true.                                      ! consider the CO2 change(need the CO2 input file)
  output_gamma = .true.                                 ! output gamma value (need gridded ef file) or real emission
  diagnose = .false.                                    ! diagnose the model of isoprene
/
```

> &#9888; Note that the variables must be adjusted to match the appropriate values for your run.

Then execute `megan_v3.3.exe`:

`> ./megan_v3.3.exe < namelist_megan` 

### Output frequency and monthly files

`output_frequency` controls the averaging interval. It defaults to `'24h'`
when omitted; values are case-insensitive. The model still runs every hour.

| Value | Output records | Records in a complete 30-day month |
| --- | --- | --- |
| `'1h'` | One hourly value per record | 720 |
| `'3h'` | Mean of each 3-hour interval | 240 |
| `'6h'` | Mean of each 6-hour interval | 120 |
| `'24h'` | Daily mean | 30 |
| `'monthly'` | Calendar-month mean | 1 |

All five choices write `Global_emis_MEGAN_YYYY-MM.nc` under `output_path`.
Records are appended within a run; a new file is created on entering each
month, including across year boundaries and leap years. Each run recreates
its monthly files, so separate runs writing the same month need separate
output directories. Static emission-factor files requested by
`output_ef_file` retain their existing annual naming.

The previous driver wrote 24 individual hourly values to each daily file.
The new default instead writes **daily means in monthly files**. Choose
`'1h'` to retain hourly resolution.

Fixed-hour intervals align with midnight UTC: for example, `'3h'` averages
00:00, 01:00 and 02:00 into the interval [00:00, 03:00). `start_date` and
`end_date` must be whole hours; `end_date` remains inclusive. To simulate
a complete month, end at 23:00 on its final day. The first or last interval
of a partial run is averaged over only the simulated hours and written
even if it is shorter than the requested frequency. A monthly mean from
a partial month therefore represents only that simulated part of the month.

The NetCDF `time` coordinate is the midpoint of each actual averaging
interval, in seconds since the first day of the file's month. `time_bnds`
contains the start and exclusive end of the interval, and `sample_count`
records the number of simulated hourly steps. Output variables have
`cell_methods = "time: mean"`. Non-finite values and `_FillValue` are excluded
at each grid cell; cells with no valid samples remain missing. `sample_count`
counts simulated steps, not the valid samples in each grid cell.

Emissions are averaged after applying the emission factors. Gamma output
and diagnostic output use the same averaging intervals. Gamma averages
are packed only after averaging, retaining the existing 0.01 scale factor;
zero is a valid gamma value, with a separate missing-value sentinel.

### Output regression tests

With a Fortran MPI compiler, NetCDF Fortran (`nf-config`), Python `numpy`
and `netCDF4` available, run from the repository root:

```sh
python -m unittest discover -s test -p 'test_output_frequency.py' -v
```

The tests compile the actual driver and run synthetic 5-longitude inputs
with one and two MPI ranks. They compare the averaged output against hourly
results, check leap-month/year boundaries and partial intervals, and cover
gamma, diagnostics, invalid options and missing-value accumulation. Override
`FC`, `FFLAGS`, `NETCDF_FFLAGS`, `NETCDF_LIBS` or `MPIEXEC` for your environment;
set `MEGAN_EXE` to use an already compiled executable. Default compiler flags
are for GNU Fortran.

Please feel free to contact the developer if you have any issues or suggestions.


---


