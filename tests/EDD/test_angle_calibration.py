from pathlib import Path
from typing import Union
import h5py
from easistrain.EDD.angleCalibEDD import angleCalibrationEDD


def generate_config(tmp_path: Path, test_data_path: Union[Path, str]) -> dict:
    HERE = Path(__file__).parent.resolve()

    with h5py.File(test_data_path, "r") as test_file:
        nb_peaks_in_boxes = test_file["infos/nbPeaksInBoxes"][()]
        fit_ranges_h = test_file["infos/rangeFitHD"][()]
        fit_ranges_v = test_file["infos/rangeFitVD"][()]
    return {
        "fileRead": str(tmp_path / "input_file.h5"),
        "fileSave": str(tmp_path / "output_file.h5"),
        "sample": "sample",
        "dataset": "0000",
        "scanNumber": 2,
        "nameHorizontalDetector": "horz_detector",
        "nameVerticalDetector": "vert_detector",
        "numberOfBoxes": len(nb_peaks_in_boxes),
        "nbPeaksInBoxes": nb_peaks_in_boxes,
        "rangeFitHD": fit_ranges_h,
        "rangeFitVD": fit_ranges_v,
        "pathFileDetectorCalibration": str(tmp_path / "calib.h5"),
        "scanDetectorCalibration": 5,
        "sampleCalibrantFile": str(HERE.parent.parent / "Calibrants" / "TiC.d"),
    }


def generate_input_files(
    tmp_path: Path,
    angle_calib_data_path: Union[Path, str],
    energy_calib_data_path: Union[Path, str],
) -> dict:
    cfg = generate_config(tmp_path, angle_calib_data_path)
    sample, dataset = cfg["sample"], cfg["dataset"]
    n_scan = cfg["scanNumber"]
    name_h, name_v = (
        cfg["nameHorizontalDetector"],
        cfg["nameVerticalDetector"],
    )
    with h5py.File(angle_calib_data_path, "r") as data_file:
        with h5py.File(cfg["fileRead"], "w") as h5file:
            h5file[f"{sample}_{dataset}_{n_scan}.1/measurement/{name_h}"] = data_file[
                "horizontal/data"
            ][()]
            h5file[f"{sample}_{dataset}_{n_scan}.1/measurement/{name_v}"] = data_file[
                "vertical/data"
            ][()]

    # Use energy calibration reference data as test input
    with h5py.File(energy_calib_data_path, "r") as data_file:
        with h5py.File(cfg["pathFileDetectorCalibration"], "w") as test_input_file:
            calib_scan_name = cfg["scanDetectorCalibration"]

            test_input_file[
                f"detectorCalibration/{calib_scan_name}/calibCoeffs/calibCoeffsHD"
            ] = data_file["horizontal/coeffs"][()]
            test_input_file[
                f"detectorCalibration/{calib_scan_name}/calibCoeffs/calibCoeffsVD"
            ] = data_file["vertical/coeffs"][()]
            test_input_file[
                f"detectorCalibration/{calib_scan_name}/calibCoeffs/uncertaintyCalibCoeffsHD"
            ] = data_file["horizontal/errors"][()]
            test_input_file[
                f"detectorCalibration/{calib_scan_name}/calibCoeffs/uncertaintyCalibCoeffsVD"
            ] = data_file["vertical/errors"][()]

    return cfg


def test_angleCalibrationEDD(tmp_path: Path):
    # Arrange
    energy_calib_data_path = (
        Path(__file__).parent.parent.resolve() / "data" / "Ba_calibration_data.hdf5"
    )
    angle_calib_data_path = (
        Path(__file__).parent.parent.resolve() / "data" / "TiC_angle_calib_data.hdf5"
    )
    with h5py.File(angle_calib_data_path, "r") as h5file:
        ref_vertical_angle = h5file["vertical/angle"][()]
        ref_horizontal_angle = h5file["horizontal/angle"][()]
        vertical_angle_error = h5file["vertical/error"][()]
        horizontal_angle_error = h5file["horizontal/error"][()]
    config = generate_input_files(
        tmp_path, angle_calib_data_path, energy_calib_data_path
    )

    # Act
    angleCalibrationEDD(**config)

    # Assert
    with h5py.File(config["fileSave"], "r") as h5file:
        grp_name = f'fit_{config["dataset"]}_{config["scanNumber"]}'
        vertical_angle = h5file[
            f"angleCalibration/{grp_name}/calibratedAngle/calibratedAngleVD"
        ][()]
        horizontal_angle = h5file[
            f"angleCalibration/{grp_name}/calibratedAngle/calibratedAngleHD"
        ][()]

    assert abs(vertical_angle - ref_vertical_angle) <= vertical_angle_error
    assert abs(horizontal_angle - ref_horizontal_angle) <= horizontal_angle_error
