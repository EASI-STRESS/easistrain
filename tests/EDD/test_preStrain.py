from pathlib import Path
import h5py
import numpy
from easistrain.EDD.preStraind0cstEDD import preStraind0cstEDD
from .utils import generate_angle_calib_file, generate_detector_calib_file


def generate_config(tmp_path: Path, test_data_path: Path) -> dict:
    with h5py.File(test_data_path, "r") as test_file:
        n_peaks = test_file["infos/nPeaks"][()]
        d0 = test_file["infos/d0"][()]

    return {
        "fileRead": str(tmp_path / "input_file.h5"),
        "fileSave": str(tmp_path / "output_file.h5"),
        "pathFileDetectorCalibration": str(tmp_path / "detector_calib.h5"),
        "scanDetectorCalibration": "scan",
        "pathFileAngleCalibration": str(tmp_path / "angle_calib.h5"),
        "scanAngleCalibration": "scan",
        "pathFileReferenceFitEDD": "",
        "scanReferenceFitEDD": "",
        "numberOfPeaks": n_peaks,
        "d0": d0,
    }


def generate_input_files(
    tmp_path: Path,
    test_data_path: Path,
    regroup_pts_data_path: Path,
    detector_data_path: Path,
    angle_data_path: Path,
) -> dict:
    cfg = generate_config(tmp_path, test_data_path)
    generate_detector_calib_file(
        detector_data_path,
        cfg["scanDetectorCalibration"],
        cfg["pathFileDetectorCalibration"],
    )
    generate_angle_calib_file(
        angle_data_path, cfg["scanAngleCalibration"], cfg["pathFileAngleCalibration"]
    )

    n_peaks = cfg["numberOfPeaks"]

    with h5py.File(cfg["fileRead"], "w") as h5file:
        with h5py.File(regroup_pts_data_path, "r") as data_file:
            for i in range(n_peaks):
                peak_grp_name = f"pointsPerPeak_{str(i).zfill(4)}"
                peak_grp = h5file.create_group(peak_grp_name)
                point_name = "point_00000"

                peak_grp[f"{point_name}"] = data_file[
                    f"peak_{i}/ref_{point_name}/value"
                ][()]
                peak_grp[f"uncertainty{point_name.capitalize()}"] = data_file[
                    f"peak_{i}/ref_{point_name}/errors"
                ][()]

    return cfg


def test_preStrain(tmp_path: Path):
    # Arrange
    data_folder = Path(__file__).parent.parent.resolve() / "data"
    test_data_path = data_folder / "BAIII_pre_strain.hdf5"
    regroup_pts_data_path = data_folder / "BAIII_regroup_points.hdf5"
    detector_data_path = data_folder / "Ba_calibration_data.hdf5"
    angle_data_path = data_folder / "TiC_angle_calib_data.hdf5"

    config = generate_input_files(
        tmp_path,
        test_data_path,
        regroup_pts_data_path,
        detector_data_path,
        angle_data_path,
    )

    # Act
    preStraind0cstEDD(**config)

    # Assert
    for i in range(config["numberOfPeaks"]):
        peak_grp_name = f"peak_{str(i).zfill(4)}"
        point_name = "point_00000"  # Take only first point
        with h5py.File(config["fileSave"], "r") as h5file:
            point_data = h5file[f"STRAIN_with_d0/{peak_grp_name}/{point_name}"][()]

        with h5py.File(test_data_path, "r") as h5file:
            ref_data = h5file[f"peak_{i}/ref_{point_name}/value"][()]
            ref_errors = h5file[f"peak_{i}/ref_{point_name}/errors"][()]

        assert numpy.all(numpy.abs(point_data - ref_data) <= numpy.abs(ref_errors))
