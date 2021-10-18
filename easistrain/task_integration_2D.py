from ewokscore import Task
from .func_integration_2D import integration_2D


class Integrate2D(
    Task,
    input_names=[
        "root_data",
        "h5file",
        "scan",
        "numScan",
        "detector_name",
        "poni_file",
        "npt_rad",
        "npt_azim",
        "x_unit",
        "im_dark",
        "im_mask",
    ],
    optional_input_names=[
        "rad_range",
        "azim_range",
        "errorModel",
        "imFlat",
        "gon1",
        "gon2",
        "gon3",
        "chiGon1",
        "omegaGon2",
        "phiGon3",
    ],
    output_names=["result"],
):
    def run(self):
        self.outputs.result = integration_2D(**self.input_values)
