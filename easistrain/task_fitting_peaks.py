from ewokscore import Task
from .func_fitting_peaks import fit


class Fit(
    Task,
    input_names=[
        "root_data",
        "h5file",
        "scan",
        "detector_name",
        "poni_file",
        "npt_rad",
        "npt_azim",
        "x_unit",
        "im_dark",
        "im_mask",
    ],
    optional_input_names=[],
    output_names=["result"],
):
    def run(self):
        self.outputs.result = fit(**self.input_values)
