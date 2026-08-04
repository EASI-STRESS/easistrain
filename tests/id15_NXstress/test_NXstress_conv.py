from easistrain.id15_NXstress import NXstressFromRaw
import os
import shutil


def test_NXstress_conv(tmp_path):
    dir_path = os.path.abspath(os.path.dirname(__file__))
    data_dir = os.path.join(dir_path, "..", "data")

    file_path = shutil.copy(
        os.path.join(data_dir, "test_id15_raw.h5"), tmp_path / "test_id15_raw.h5"
    )
    det_calib_file_angle = shutil.copy(
        os.path.join(data_dir, "angleCalib.h5"), tmp_path / "angleCalib.h5"
    )
    det_calib_file_energy = shutil.copy(
        os.path.join(data_dir, "energyCalib.h5"), tmp_path / "energyCalib.h5"
    )

    nx_stress = NXstressFromRaw(
        file_path=file_path,
        det_calib_file_angle=det_calib_file_angle,
        det_calib_file_energy=det_calib_file_energy,
        with_cradle=False,
        lattice="lattice",
        phase_name="phase_name",
        scanNbForRotation=2,
        experimental_identifier="experimental_identifier",
        collection_identifier="collection_identifier",
        test_script=True,
    )
    nx_stress.main()
