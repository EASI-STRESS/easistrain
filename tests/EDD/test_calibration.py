from pathlib import Path
from typing import Union
import numpy
import h5py
from easistrain.EDD.calibrationEDD import calibEdd as calib_edd


def test_calib_edd(tmp_path: Path):
    test_data_path, config = calib_edd_init(tmp_path)
    calib_edd(**config)
    calib_edd_assert(test_data_path, config)


def generate_config(tmp_path: Path, test_data_path: Union[Path, str]) -> dict:
    HERE = Path(__file__).parent.resolve()

    with h5py.File(test_data_path, "r") as test_file:
        nb_peaks_in_boxes = test_file["infos/nbPeaksInBoxes"][()]
        fit_ranges = test_file["infos/rangeFit"][()]
    return {
        "fileRead": str(tmp_path / "input_file.h5"),
        "fileSave": str(tmp_path / "output_file.h5"),
        "sample": "sample",
        "dataset": "0000",
        "scanNumberHorizontalDetector": 2,
        "scanNumberVerticalDetector": 1,
        "nameHorizontalDetector": "horz_detector",
        "nameVerticalDetector": "vert_detector",
        "nbPeaksInBoxes": nb_peaks_in_boxes,
        "rangeFit": fit_ranges,
        "sourceCalibrantFile": str(HERE.parent.parent / "Calibrants" / "BaSource"),
    }


def generate_input_files(tmp_path: Path, test_data_path: Union[Path, str]) -> dict:
    cfg = generate_config(tmp_path, test_data_path)

    sample, dataset = cfg["sample"], cfg["dataset"]
    n_scan_h, name_h = (
        cfg["scanNumberHorizontalDetector"],
        cfg["nameHorizontalDetector"],
    )
    n_scan_v, name_v = (
        cfg["scanNumberVerticalDetector"],
        cfg["nameVerticalDetector"],
    )
    with h5py.File(test_data_path, "r") as test_file:
        with h5py.File(cfg["fileRead"], "w") as h5file:
            h5file[f"{sample}_{dataset}_{n_scan_h}.1/measurement/{name_h}"] = test_file[
                "horizontal/data"
            ][()]
            h5file[f"{sample}_{dataset}_{n_scan_v}.1/measurement/{name_v}"] = test_file[
                "vertical/data"
            ][()]

    return cfg


def calib_edd_init(tmp_path: Path):
    test_data_path = (
        Path(__file__).parent.parent.resolve() / "data" / "Ba_calibration_data.hdf5"
    )
    config = generate_input_files(tmp_path, test_data_path)
    return test_data_path, config


def calib_edd_assert(test_data_path, config):
    with h5py.File(config["fileSave"], "r") as h5file:
        grp_name = f'fit_{config["dataset"]}_{config["scanNumberHorizontalDetector"]}_{config["scanNumberVerticalDetector"]}'
        vertical_coeffs = h5file[
            f"detectorCalibration/{grp_name}/calibCoeffs/calibCoeffsVD"
        ][()]
        horizontal_coeffs = h5file[
            f"detectorCalibration/{grp_name}/calibCoeffs/calibCoeffsHD"
        ][()]

    with h5py.File(test_data_path, "r") as h5file:
        ref_vertical_coeffs = h5file["vertical/coeffs"][()]
        ref_horizontal_coeffs = h5file["horizontal/coeffs"][()]
        vertical_coeffs_errors = h5file["vertical/errors"][()]
        horizontal_coeffs_errors = h5file["horizontal/errors"][()]

    assert numpy.all(
        numpy.abs(vertical_coeffs - ref_vertical_coeffs) <= vertical_coeffs_errors
    )
    assert numpy.all(
        numpy.abs(horizontal_coeffs - ref_horizontal_coeffs) <= horizontal_coeffs_errors
    )
