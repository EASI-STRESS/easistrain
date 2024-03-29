from pathlib import Path
from typing import Optional, Union
import h5py
import numpy
import pytest
from easistrain.EDD.fitEDD import fitEDD

ORIENTATION = "OR1"


def generate_config(
    tmp_path: Path, test_data_path: Union[Path, str], scan_length: Optional[int] = None
) -> dict:
    with h5py.File(test_data_path, "r") as test_file:
        nb_peaks_in_boxes = test_file[f"{ORIENTATION}/infos/nbPeaksInBoxes"][()]
        fit_ranges_h = test_file[f"{ORIENTATION}/infos/rangeFitHD"][()]
        fit_ranges_v = test_file[f"{ORIENTATION}/infos/rangeFitVD"][()]
        positioners = list(test_file[f"{ORIENTATION}/positioners"].keys())

    return {
        "fileRead": str(tmp_path / "input_file.h5"),
        "fileSave": str(tmp_path / "output_file.h5"),
        "sample": "sample",
        "dataset": "0000",
        "scanNumber": 2,
        "nameHorizontalDetector": "horz_detector",
        "nameVerticalDetector": "vert_detector",
        "positioners": positioners,
        "nbPeaksInBoxes": nb_peaks_in_boxes,
        "rangeFitHD": fit_ranges_h,
        "rangeFitVD": fit_ranges_v,
        "detectorSliceIndex": None if scan_length is None else scan_length // 2,
    }


def generate_input_files(
    tmp_path: Path, test_data_path: Union[Path, str], scan_length: Optional[int] = None
) -> dict:
    cfg = generate_config(tmp_path, test_data_path, scan_length)
    sample, dataset = cfg["sample"], cfg["dataset"]
    n_scan = cfg["scanNumber"]
    name_h, name_v = (
        cfg["nameHorizontalDetector"],
        cfg["nameVerticalDetector"],
    )

    if scan_length:
        detector_dim = scan_length
    else:
        detector_dim = 1

    with h5py.File(test_data_path, "r") as data_file:
        with h5py.File(cfg["fileRead"], "w") as h5file:
            scan_grp = f"{sample}_{dataset}_{n_scan}.1"
            raw_horz_dataset = data_file[f"{ORIENTATION}/horizontal/data"]
            assert isinstance(raw_horz_dataset, h5py.Dataset)
            horz_dataset = numpy.zeros(
                (detector_dim, len(raw_horz_dataset)), dtype=raw_horz_dataset.dtype
            )
            for i in range(horz_dataset.shape[0]):
                horz_dataset[i] = raw_horz_dataset[()]
            h5file.create_dataset(
                f"{scan_grp}/measurement/{name_h}",
                data=horz_dataset,
            )
            raw_vert_dataset = data_file[f"{ORIENTATION}/vertical/data"]
            assert isinstance(raw_vert_dataset, h5py.Dataset)
            vert_dataset = numpy.zeros(
                (detector_dim, len(raw_vert_dataset)), dtype=raw_vert_dataset.dtype
            )
            for i in range(vert_dataset.shape[0]):
                vert_dataset[i] = raw_vert_dataset[()]
            h5file.create_dataset(
                f"{scan_grp}/measurement/{name_v}",
                data=vert_dataset,
            )
            for pos_name, pos_data in data_file[f"{ORIENTATION}/positioners"].items():
                h5file[f"{scan_grp}/instrument/positioners/{pos_name}"] = pos_data[()]

    return cfg


@pytest.mark.parametrize("scan_length", [None, 5])
def test_fitEDD(tmp_path: Path, scan_length: Optional[int]):
    # Arrange
    test_data_path = (
        Path(__file__).parent.parent.resolve() / "data" / "BAIII_fit_data.hdf5"
    )
    with h5py.File(test_data_path, "r") as h5file:
        ref_vertical_params = h5file[f"{ORIENTATION}/vertical/fit/params"][()]
        vertical_params_errors = h5file[f"{ORIENTATION}/vertical/fit/errors"][()]
        ref_horizontal_params = h5file[f"{ORIENTATION}/horizontal/fit/params"][()]
        horizontal_params_errors = h5file[f"{ORIENTATION}/horizontal/fit/errors"][()]
    config = generate_input_files(tmp_path, test_data_path, scan_length)

    # Act
    fitEDD(**config)

    # Assert
    with h5py.File(config["fileSave"], "r") as h5file:
        grp_name = (
            f'{config["sample"]}_{config["dataset"]}_{config["scanNumber"]}.1/fit'
        )
        vertical_params = h5file[f"{grp_name}/0000/fitParams/fitParamsVD"][()]
        horizontal_params = h5file[f"{grp_name}/0000/fitParams/fitParamsHD"][()]

    # Remove latest column that is the goodness factor
    assert numpy.all(
        numpy.abs(vertical_params - ref_vertical_params)[:, :-1]
        <= vertical_params_errors
    )
    assert numpy.all(
        numpy.abs(horizontal_params - ref_horizontal_params)[:, :-1]
        <= horizontal_params_errors
    )
