from easistrain.id15_NXstress import NXstressFromRaw
import os


def test_NXstress_conv():
    dir_path = os.path.abspath(os.path.dirname(__file__))

    file_path = os.path.join(dir_path, "..", "data", "test_id15_raw.h5")
    det_calib_file_angle = os.path.join(dir_path, "..", "data", "angleCalib.h5")
    det_calib_file_energy = os.path.join(dir_path, "..", "data", "energyCalib.h5")

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
