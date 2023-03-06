from typing import Callable, Dict, Sequence, Union
from typing_extensions import Literal

import h5py
import numpy

from easistrain.EDD.io import save_fit_data, save_fit_params
from easistrain.EDD.utils import fit_detector_data

Detector = Literal["horizontal", "vertical"]

DETECTORS: Sequence[Detector] = ["horizontal", "vertical"]


def fit_all_peaks_and_save_results(
    nbPeaksInBoxes: Sequence[int],
    rangeFit: Dict[Detector, Sequence[int]],
    patterns: Dict[Detector, numpy.ndarray],
    scanNumbers: Dict[Detector, Union[str, int]],
    saving_dest: h5py.Group,
    group_format: Callable[[int], str],
):
    fitParams = {"horizontal": numpy.array(()), "vertical": numpy.array(())}
    uncertaintyFitParams = {
        "horizontal": numpy.array(()),
        "vertical": numpy.array(()),
    }
    for i, nb_peaks in enumerate(nbPeaksInBoxes):
        fit_line_group = saving_dest.create_group(
            group_format(i)
        )  ## create group for each calibration peak

        for detector in DETECTORS:
            fit_min, fit_max = (
                rangeFit[detector][2 * i],
                rangeFit[detector][2 * i + 1],
            )
            scanNumber = scanNumbers[detector]
            channels = numpy.arange(fit_min, fit_max)
            raw_data = patterns[detector][fit_min:fit_max]
            assert isinstance(raw_data, numpy.ndarray)

            (
                background,
                fitted_data,
                boxFitParams,
                uncertaintyBoxFitParams,
            ) = fit_detector_data(
                channels=channels,
                raw_data=raw_data,
                nb_peaks=nb_peaks,
                boxCounter=i,
                scanNumber=int(scanNumber),
                detectorName=detector,
            )

            save_fit_data(
                fit_line_group, detector, channels, raw_data, background, fitted_data
            )

            # Accumulate fit parameters of this box
            fitParams[detector] = numpy.append(fitParams[detector], boxFitParams)
            uncertaintyFitParams[detector] = numpy.append(
                uncertaintyFitParams[detector], uncertaintyBoxFitParams
            )

    fit_params_group = saving_dest.create_group("fitParams")
    return save_fit_params(fit_params_group, fitParams, uncertaintyFitParams)
