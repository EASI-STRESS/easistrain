from pathlib import Path
import h5py
import numpy
from easistrain.EDD.strainStressd0cstEDD import strainStressTensor


def generate_config(tmp_path: Path, test_data_path: Path) -> dict:
    with h5py.File(test_data_path, "r") as test_file:
        n_peaks = test_file["infos/nPeaks"][()]
        XEC = test_file["infos/XEC"][()]

    return {
        "fileRead": str(tmp_path / "input_file.h5"),
        "fileSave": str(tmp_path / "output_file.h5"),
        "numberOfPeaks": n_peaks,
        "XEC": XEC,
    }


def generate_input_files(
    tmp_path: Path,
    test_data_path: Path,
    pre_strain_data_path: Path,
) -> dict:
    cfg = generate_config(tmp_path, test_data_path)
    n_peaks = cfg["numberOfPeaks"]

    with h5py.File(cfg["fileRead"], "w") as h5file:
        with h5py.File(pre_strain_data_path, "r") as data_file:
            point_name = "point_00000"
            for i in range(n_peaks):
                peak_grp = h5file.create_group(f"STRAIN_with_d0/peak_{str(i).zfill(4)}")

                peak_grp[f"{point_name}"] = data_file[
                    f"peak_{i}/ref_{point_name}/value"
                ][()]
                peak_grp[f"uncertainty_{point_name}"] = data_file[
                    f"peak_{i}/ref_{point_name}/errors"
                ][()]

    return cfg


def test_strain(tmp_path: Path):
    # Arrange
    data_folder = Path(__file__).parent.parent.resolve() / "data"
    test_data_path = data_folder / "BAIII_strain.hdf5"
    pre_strain_data_path = data_folder / "BAIII_pre_strain.hdf5"

    config = generate_input_files(
        tmp_path,
        test_data_path,
        pre_strain_data_path,
    )

    # Act
    strainStressTensor(**config)

    # Assert
    for tensor_type in ["strain", "stress"]:
        for i in range(config["numberOfPeaks"]):
            peak_grp_name = f"peak_{str(i).zfill(4)}"
            point_name = "point_00000"  # Take only first point
            with h5py.File(config["fileSave"], "r") as h5file:
                point_data = h5file[
                    f"{tensor_type}_tensor/{peak_grp_name}/{point_name}"
                ][()]

            with h5py.File(test_data_path, "r") as h5file:
                ref_data = h5file[f"{tensor_type}/peak_{i}/ref_{point_name}/value"][()]
                ref_errors = h5file[f"{tensor_type}/peak_{i}/ref_{point_name}/errors"][
                    ()
                ]

            assert numpy.all(numpy.abs(point_data - ref_data) <= numpy.abs(ref_errors))
