from typing import Union
from pathlib import Path
import h5py


def generate_detector_calib_file(
    test_data_path: Union[Path, str], scan_name: str, output_path: Union[Path, str]
):
    with h5py.File(test_data_path, "r") as test_file:
        with h5py.File(output_path, "w") as h5file:
            h5file[
                f"detectorCalibration/{scan_name}/calibCoeffs/calibCoeffsVD"
            ] = test_file["vertical/coeffs"][()]
            h5file[
                f"detectorCalibration/{scan_name}/calibCoeffs/uncertaintyCalibCoeffsVD"
            ] = test_file["vertical/errors"][()]

            h5file[
                f"detectorCalibration/{scan_name}/calibCoeffs/calibCoeffsHD"
            ] = test_file["horizontal/coeffs"][()]

            h5file[
                f"detectorCalibration/{scan_name}/calibCoeffs/uncertaintyCalibCoeffsHD"
            ] = test_file["horizontal/errors"][()]


def generate_angle_calib_file(
    test_data_path: Union[Path, str], scan_name: str, output_path: Union[Path, str]
):
    with h5py.File(test_data_path, "r") as test_file:
        with h5py.File(output_path, "w") as h5file:
            h5file[
                f"angleCalibration/{scan_name}/calibratedAngle/calibratedAngleVD"
            ] = test_file["vertical/angle"][()]
            h5file[
                f"angleCalibration/{scan_name}/calibratedAngle/uncertaintyCalibratedAngleVD"
            ] = test_file["vertical/error"][()]

            h5file[
                f"angleCalibration/{scan_name}/calibratedAngle/calibratedAngleHD"
            ] = test_file["horizontal/angle"][()]

            h5file[
                f"angleCalibration/{scan_name}/calibratedAngle/uncertaintyCalibratedAngleHD"
            ] = test_file["horizontal/error"][()]
