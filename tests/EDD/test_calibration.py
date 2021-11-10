from pathlib import Path
from typing import Union
import numpy as np
import h5py
from easistrain.EDD.calibrationEDD import calibEdd


def generate_config(tmp_path: Path) -> dict:
    HERE = Path(__file__).parent.resolve()
    return {
        "fileRead": str(tmp_path / "input_file.h5"),
        "fileSave": str(tmp_path / "output_file.h5"),
        "sample": "sample",
        "dataset": "0000",
        "scanNumberHorizontalDetector": 2,
        "scanNumberVerticalDetector": 1,
        "nameHorizontalDetector": "horz_detector",
        "nameVerticalDetector": "vert_detector",
        "numberOfBoxes": 4,
        "nbPeaksInBoxes": [1, 2, 1, 1],
        "rangeFit": [620, 780, 1020, 1120, 3500, 3800, 3850, 4090],
        "sourceCalibrantFile": str(HERE.parent.parent / "Calibrants" / "BaSource"),
    }


def generate_input_file(cfg: dict, test_data_path: Union[Path, str]):
    sample, dataset = cfg["sample"], cfg["dataset"]
    n_scan_h, name_h = (
        cfg["scanNumberHorizontalDetector"],
        cfg["nameHorizontalDetector"],
    )
    n_scan_v, name_v = (
        cfg["scanNumberVerticalDetector"],
        cfg["nameVerticalDetector"],
    )
    with h5py.File(test_data_path, "r") as testFile:
        with h5py.File(cfg["fileRead"], "w") as h5file:
            h5file[f"{sample}_{dataset}_{n_scan_h}.1/measurement/{name_h}"] = testFile[
                "horizontal/data"
            ][()]
            h5file[f"{sample}_{dataset}_{n_scan_v}.1/measurement/{name_v}"] = testFile[
                "vertical/data"
            ][()]


def test_calibEdd(tmp_path: Path):
    # Arrange
    config = generate_config(tmp_path)
    test_data_path = (
        Path(__file__).parent.parent.resolve() / "data" / "Ba_calibration_data.hdf5"
    )
    with h5py.File(test_data_path, "r") as h5file:
        ref_vertical_coeffs = h5file["vertical/coeffs"][()]
        ref_horizontal_coeffs = h5file["horizontal/coeffs"][()]
        vertical_coeffs_errors = h5file["vertical/errors"][()]
        horizontal_coeffs_errors = h5file["horizontal/errors"][()]
    generate_input_file(config, test_data_path)

    # Act
    calibEdd(**config)

    # Assert
    with h5py.File(config["fileSave"], "r") as h5file:
        grp_name = f'fit_{config["dataset"]}_{config["scanNumberHorizontalDetector"]}_{config["scanNumberVerticalDetector"]}'
        vertical_coeffs = h5file[
            f"detectorCalibration/{grp_name}/calibCoeffs/calibCoeffsVD"
        ][()]
        horizontal_coeffs = h5file[
            f"detectorCalibration/{grp_name}/calibCoeffs/calibCoeffsHD"
        ][()]

    assert np.all(
        np.abs(vertical_coeffs - ref_vertical_coeffs) <= vertical_coeffs_errors
    )
    assert np.all(
        np.abs(horizontal_coeffs - ref_horizontal_coeffs) <= horizontal_coeffs_errors
    )
