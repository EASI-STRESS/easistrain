from pathlib import Path
from typing import Union
import h5py
import numpy
import os.path
from easistrain.EDD.regroupPoints import regroupPoints

ORIENTATIONS = ["OR1", "OR2"]


def generate_config(tmp_path: Path, test_data_path: Union[Path, str]) -> dict:
    with h5py.File(test_data_path, "r") as test_file:
        n_peaks = test_file["infos/nPeaks"][()]

    return {
        "fileRead": [
            str(tmp_path / f"{orientation}.h5") for orientation in ORIENTATIONS
        ],
        "fileSave": str(tmp_path / "output_file.h5"),
        "numberOfPeaks": n_peaks,
    }


def generate_input_files(
    tmp_path: Path,
    test_data_path: Union[Path, str],
    coord_transform_data_path: Union[Path, str],
) -> dict:
    cfg = generate_config(tmp_path, test_data_path)
    n_peaks = cfg["numberOfPeaks"]

    for orientation_file in cfg["fileRead"]:
        orientation = os.path.basename(orientation_file).split(".")[0]
        with h5py.File(orientation_file, "w") as h5file:
            global_grp = h5file.create_group("global")
            with h5py.File(coord_transform_data_path, "r") as data_file:
                for i in range(n_peaks):
                    peak_name = f"peak_{str(i).zfill(4)}"
                    global_grp[f"inSample_{peak_name}"] = data_file[
                        f"{orientation}/ref_{peak_name}/value"
                    ][()]
                    global_grp[
                        f"inSample_uncertainty{peak_name.capitalize()}"
                    ] = data_file[f"{orientation}/ref_{peak_name}/errors"][()]

    return cfg


def test_regroupPoints(tmp_path: Path):
    # Arrange
    data_folder = Path(__file__).parent.parent.resolve() / "data"
    test_data_path = data_folder / "BAIII_regroup_points.hdf5"
    coord_transform_data_path = data_folder / "BAIII_coord_transform.hdf5"

    config = generate_input_files(tmp_path, test_data_path, coord_transform_data_path)

    # Act
    regroupPoints(**config)

    # Assert
    for i in range(config["numberOfPeaks"]):
        peak_grp_name = f"pointsPerPeak_{str(i).zfill(4)}"
        point_name = "point_00000"  # Take only first point
        with h5py.File(config["fileSave"], "r") as h5file:
            point_data = h5file[f"{peak_grp_name}/{point_name}"][()]

        with h5py.File(test_data_path, "r") as h5file:
            ref_data = h5file[f"peak_{i}/ref_{point_name}/value"][()]
            ref_errors = h5file[f"peak_{i}/ref_{point_name}/errors"][()]

        print(numpy.abs(point_data - ref_data))
        assert numpy.all(numpy.abs(point_data - ref_data) <= numpy.abs(ref_errors))
