"""Integration tests using tiny synthetic inputs and the real MPI/NetCDF driver.

Requires numpy, netCDF4, a Fortran MPI compiler, and NetCDF Fortran.
Run: python -m unittest discover -s test -p 'test_output_frequency.py' -v
Optional overrides: MEGAN_EXE, FC, FFLAGS, NETCDF_FFLAGS, NETCDF_LIBS, MPIEXEC.
The default build flags target GNU Fortran; set FFLAGS for other compilers.
"""

import calendar
from datetime import datetime, timedelta
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest

import netCDF4
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
HOUR = timedelta(hours=1)


def flags(variable, default):
    if variable in os.environ:
        return shlex.split(os.environ[variable])
    return shlex.split(subprocess.check_output(default, text=True))


def run_checked(command, **kwargs):
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, timeout=180, **kwargs)
    if result.returncode:
        raise AssertionError(f"{command}\n{result.stdout[-8000:]}")
    return result


def grid_file(path):
    ds = netCDF4.Dataset(path, "w")
    ds.createDimension("latitude", 5)
    ds.createDimension("longitude", 5)
    ds.createVariable("latitude", "f4", ("latitude",))[:] = [90, 30, 0, -30, -90]
    ds.createVariable("longitude", "f4", ("longitude",))[:] = [0, 15, 30, 45, 60]
    return ds


def make_inputs(path):
    """Five longitudes exercise uneven decomposition across two MPI ranks."""
    xy = np.arange(25, dtype=np.float32).reshape(5, 5) / 10
    for year, month in [(2023, 11), (2023, 12), (2024, 1), (2024, 2), (2024, 3)]:
        start = datetime(year, month, 1)
        hours = calendar.monthrange(year, month)[1] * 24
        hour = np.arange(hours)[:, None, None] % 24
        with grid_file(path / f"met_{year}_{month:02d}.nc") as ds:
            ds.createDimension("valid_time", hours)
            ds.createVariable("valid_time", "i8", ("valid_time",))[:] = [
                int((start + i * HOUR - datetime(1970, 1, 1)).total_seconds()) for i in range(hours)
            ]
            fields = {
                "t2m": 294 + 4 * np.sin(hour * np.pi / 12) + xy,
                "d2m": 283 + np.zeros_like(hour) + xy,
                "u10": 2 + np.zeros_like(hour) + xy / 5,
                "v10": 1 + np.zeros_like(hour) + xy / 10,
                "sp": 101325 + np.zeros_like(hour) + xy,
                # ERA5-Land accumulated radiation resets at 01:00.
                "ssrd": np.where(hour == 0, 24, hour) * 3600 * (180 + xy),
            }
            for name, data in fields.items():
                ds.createVariable(name, "f4", ("valid_time", "latitude", "longitude"))[:] = data
    for year in (2023, 2024):
        with grid_file(path / f"pft_{year}.nc") as ds:
            for name in ("NEEDTR", "TROPBE", "BROADTR", "SHRUB", "GRASS", "CROP", "TREE"):
                data = np.zeros((5, 5), dtype=np.float32)
                if name == "SHRUB":
                    data[:, 1:] = 1.0  # Longitude zero has no vegetation: valid zero emissions/gamma.
                ds.createVariable(name, "f4", ("latitude", "longitude"))[:] = data
    for year in (2022, 2023, 2024):
        for index in range(1, 13):
            with grid_file(path / f"lai_{year}{index:03d}.nc") as ds:
                ds.createVariable("LAI", "i1", ("latitude", "longitude"))[:] = 30 + index
    with grid_file(path / "ef.nc") as ds:
        ds.createDimension("species", 19)
        ds.createDimension("ldf_species", 4)
        for group in ("SHRUB", "GRASS", "CROP", "TREE"):
            ds.createVariable(group + "_EF", "f4", ("species", "latitude", "longitude"))[:] = (
                np.arange(1, 20)[:, None, None] * (1 + xy)
            )
            ds.createVariable(group + "_LDF", "f4", ("ldf_species", "latitude", "longitude"))[:] = 0.6


class OutputFrequencyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix="megan-output-test-")
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.path = Path(cls.tmp.name)
        cls.inputs = cls.path / "inputs"
        cls.inputs.mkdir()
        make_inputs(cls.inputs)
        cls.compiler = shlex.split(os.environ.get("FC", "mpifort"))
        cls.fflags = shlex.split(os.environ.get("FFLAGS", "-O0 -g -ffree-line-length-none -fcheck=bounds"))
        cls.ncflags = flags("NETCDF_FFLAGS", ["nf-config", "--fflags"])
        cls.nclibs = flags("NETCDF_LIBS", ["nf-config", "--flibs"])
        cls.launch = shlex.split(os.environ.get("MPIEXEC", "mpiexec"))
        cls.exe = Path(os.environ.get("MEGAN_EXE", cls.path / "megan.exe")).resolve()
        if "MEGAN_EXE" not in os.environ:
            sources = [ROOT / "src" / name for name in (
                "datetime_module.f90", "megan_voc.f90", "read_co2.f90", "megan_driver.f90")]
            run_checked(cls.compiler + cls.fflags + cls.ncflags + ["-I" + str(ROOT / "src")]
                        + list(map(str, sources)) + cls.nclibs + ["-o", str(cls.exe)], cwd=cls.path)
        cls.run_number = 0

    def simulate(self, frequency, start, end, ranks=2, gamma=False, diagnose=False, success=True):
        type(self).run_number += 1
        out = self.path / f"run_{self.run_number}"
        out.mkdir()
        frequency_line = "" if frequency is None else f"output_frequency='{frequency}'"
        namelist = f"""&megan_nl
start_date='{start}', end_date='{end}'
{frequency_line}
met_file_path='{self.inputs}', met_files='met_<time>.nc'
pft_file_path='{self.inputs}', pft_files='pft_<time>.nc'
lai_file_path='{self.inputs}', lai_files='lai_<time>.nc'
ef_file_path='{self.inputs}', ef_file='ef.nc'
output_path='{out}'
nlai='12', lai_scale_factor=0.1
run_co2=.false., real_spinup_met=.false.
output_gamma={'.true.' if gamma else '.false.'}
diagnose={'.true.' if diagnose else '.false.'}
/
"""
        command = self.launch + ["-n", str(ranks), str(self.exe)]
        result = subprocess.run(command, input=namelist, text=True, cwd=self.path,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        if not success:
            self.assertNotEqual(result.returncode, 0, result.stdout[-4000:])
            self.assertFalse(list(out.glob("*.nc")))
            return result.stdout
        self.assertEqual(result.returncode, 0, result.stdout[-6000:])
        records = {}
        for file in sorted(out.glob("*.nc")):
            self.assertRegex(file.name, r"^Global_emis_MEGAN_\d{4}-\d{2}\.nc$")
            with netCDF4.Dataset(file) as ds:
                self.assertTrue(ds.dimensions["time"].isunlimited())
                self.assertEqual(ds.output_frequency, (frequency or "24h").strip().lower())
                time = ds.variables["time"]
                self.assertEqual(time.bounds, "time_bnds")
                origin = datetime.strptime(file.stem[-7:] + "-01", "%Y-%m-%d")
                bounds = ds.variables["time_bnds"][:]
                np.testing.assert_array_equal(time[:], bounds.mean(axis=1))
                # Ensure the declared units decode to the same actual midpoint.
                dates = netCDF4.num2date(time[:], time.units, time.calendar)
                for i, (left, right) in enumerate(bounds):
                    first, next_ = origin + timedelta(seconds=float(left)), origin + timedelta(seconds=float(right))
                    self.assertEqual(dates[i].strftime("%Y-%m-%d %H:%M:%S"),
                                     (first + (next_ - first) / 2).strftime("%Y-%m-%d %H:%M:%S"))
                    self.assertEqual(int(ds.variables["sample_count"][i]), int((next_ - first) / HOUR))
                    fields = {}
                    for name, var in ds.variables.items():
                        if var.dimensions == ("time", "lat", "lon"):
                            self.assertEqual(var.cell_methods, "time: mean")
                            data = var[i]
                            # Zero-vegetation cells must be valid zeros, including packed gamma.
                            self.assertFalse(np.ma.getmaskarray(data)[:, 0].any())
                            np.testing.assert_array_equal(data[:, 0], 0)
                            fields[name] = data.filled(np.nan)
                    self.assertEqual(len(fields), 19)
                    self.assertNotIn((first, next_), records)
                    records[first, next_] = fields
        self.assertTrue(records)
        return records

    def compare_means(self, hourly, averaged, tolerance=0):
        consumed = 0
        for (start, end), fields in averaged.items():
            subset = [values for (a, b), values in hourly.items() if start <= a and b <= end]
            self.assertEqual(len(subset), int((end - start) / HOUR))
            consumed += len(subset)
            for name, actual in fields.items():
                expected = np.mean([values[name] for values in subset], axis=0, dtype=np.float64)
                np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=tolerance, equal_nan=True)
        self.assertEqual(consumed, len(hourly))

    def test_all_frequencies_partial_windows_and_leap_day(self):
        start, end = "2024-02-28 01:00:00", "2024-03-01 04:00:00"
        hourly = self.simulate("1h", start, end, ranks=1)
        self.assertEqual(len(hourly), 52)
        for frequency, count in [("3h", 18), ("6h", 9), ("24h", 3), ("monthly", 2)]:
            with self.subTest(frequency=frequency):
                averaged = self.simulate(frequency, start, end)
                self.assertEqual(len(averaged), count)
                self.compare_means(hourly, averaged)
                for (_, right) in list(averaged)[:-1]:
                    if frequency != "monthly":
                        self.assertEqual(right.hour % int(frequency[:-1]), 0)

    def test_complete_leap_month(self):
        start, end = "2024-02-01 00:00:00", "2024-02-29 23:00:00"
        daily = self.simulate("24h", start, end, ranks=1)
        monthly = self.simulate("monthly", start, end)
        self.assertEqual(len(daily), 29)
        self.assertEqual(list(monthly), [(datetime(2024, 2, 1), datetime(2024, 3, 1))])
        for name, actual in next(iter(monthly.values())).items():
            np.testing.assert_allclose(actual, np.mean([d[name] for d in daily.values()], axis=0,
                                                       dtype=np.float64), rtol=2e-6, equal_nan=True)

    def test_year_rollover(self):
        start, end = "2023-12-31 20:00:00", "2024-01-01 04:00:00"
        hourly = self.simulate("1h", start, end)
        monthly = self.simulate("MONTHLY", start, end)
        self.assertEqual(len(monthly), 2)
        self.compare_means(hourly, monthly)

    def test_gamma_and_diagnostics(self):
        for mode in ({"gamma": True}, {"diagnose": True, "gamma": True}):
            with self.subTest(mode=mode):
                start, end = "2024-02-29 01:00:00", "2024-03-01 04:00:00"
                hourly = self.simulate("1h", start, end, ranks=1, **mode)
                averaged = self.simulate("6H", start, end, **mode)
                # Packed hourly gamma differs from packing after averaging by < 0.01.
                self.compare_means(hourly, averaged, tolerance=0 if mode.get("diagnose") else 0.010001)

    def test_default_single_hour_and_invalid_settings(self):
        start, end = "2024-02-29 01:00:00", "2024-02-29 01:00:00"
        hourly = self.simulate("1h", start, end)
        for frequency in (None, "monthly", "3h", "6h", "24h"):
            self.compare_means(hourly, self.simulate(frequency, start, end))
        for frequency in ("2h", "daily", ""):
            message = self.simulate(frequency, start, end, success=False)
            self.assertIn("output_frequency must be", message)
        self.assertIn("precedes", self.simulate("1h", end, "2024-02-28 01:00:00", success=False))
        self.assertIn("whole hours", self.simulate("1h", "2024-02-29 01:30:00",
                                                   "2024-02-29 02:00:00", success=False))

    def test_missing_values_in_accumulator(self):
        # Compile the actual helper, including NaN and fill inputs that the tiny
        # vegetation fixture does not otherwise generate.
        driver = (ROOT / "src/megan_driver.f90").read_text()
        helper = re.search(r"(?ms)^   subroutine accumulate_output\(.*?^   end subroutine accumulate_output", driver)[0]
        program = """program test_accumulator
use iso_fortran_env, only: real64
use ieee_arithmetic
implicit none
real, parameter :: FillValue=-999.
real :: values(4,1,1)
real(real64) :: sums(4,1,1)=0
integer :: counts(4,1,1)=0
values(:,1,1)=[1., FillValue, ieee_value(0.,ieee_quiet_nan), 0.]
call accumulate_output(values,sums,counts)
values(:,1,1)=[3., 5., FillValue, 0.]
call accumulate_output(values,sums,counts)
if(any(counts(:,1,1) /= [2,1,0,2])) stop 1
if(any(sums(:,1,1) /= [4._real64,5._real64,0._real64,0._real64])) stop 2
contains
""" + helper + "\nend program\n"
        source = self.path / "test_accumulator.f90"
        source.write_text(program)
        exe = self.path / "test_accumulator"
        run_checked(self.compiler + self.fflags + [str(source), "-o", str(exe)], cwd=self.path)
        run_checked([str(exe)], cwd=self.path)


if __name__ == "__main__":
    unittest.main(verbosity=2)
