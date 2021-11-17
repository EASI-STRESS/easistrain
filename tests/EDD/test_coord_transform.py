from pathlib import Path
from typing import Union
import h5py
import numpy
from easistrain.EDD.coordTransformation import coordTransformation


def generate_config(tmp_path: Path, test_data_path: Union[Path, str]) -> dict:
    with h5py.File(test_data_path, "r") as test_file:
        n_peaks = test_file["infos/nPeaks"][()]
        gonioToSample = test_file["infos/gonioToSample"][()]

    return {
        "fileRead": str(tmp_path / "input_file.h5"),
        "fileSave": str(tmp_path / "output_file.h5"),
        "numberOfPeaks": n_peaks,
        "gonioToSample": gonioToSample.tolist(),
    }


def generate_input_files(
    tmp_path: Path,
    test_data_path: Union[Path, str],
) -> dict:
    cfg = generate_config(tmp_path, test_data_path)
    n_peaks = cfg["numberOfPeaks"]

    with h5py.File(test_data_path, "r") as data_file:
        with h5py.File(cfg["fileRead"], "w") as h5file:
            pos_group = h5file.create_group("scan/tthPositionsGroup")
            for i in range(n_peaks):
                peak_name = f"peak_{str(i).zfill(4)}"
                print(test_data_path, list(data_file.keys()))
                pos_group[peak_name] = data_file[f"{peak_name}/value"][()]
                pos_group[f"uncertainty{peak_name.capitalize()}"] = data_file[
                    f"{peak_name}/errors"
                ][()]

    return cfg


def test_coordTransform(tmp_path: Path):
    # Arrange
    test_data_path = (
        Path(__file__).parent.parent.resolve() / "data" / "BAIII_coord_transform.hdf5"
    )
    config = generate_input_files(tmp_path, test_data_path)

    # Act
    coordTransformation(**config)

    # Assert
    for i in range(config["numberOfPeaks"]):
        peak_name = f"peak_{str(i).zfill(4)}"
        with h5py.File(config["fileSave"], "r") as h5file:
            peak_data = h5file[f"global/inSample_{peak_name}"][()]

        with h5py.File(test_data_path, "r") as h5file:
            ref_data = h5file[f"ref_{peak_name}/value"][()]
            ref_errors = h5file[f"ref_{peak_name}/errors"][()]

        assert numpy.all(numpy.abs(peak_data - ref_data) <= numpy.abs(ref_errors))
