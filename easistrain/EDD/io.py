import numpy as np


def peak_dataset_data(positionAngles: np.ndarray, savedPeakFitParams: np.ndarray):
    return [
        *positionAngles[0, 0:6],
        -90,  ## delta angle of the horizontal detector (debye scherer angle)
        0,  ## theta angle (diffraction fixed angle) of the horizontal detector (I suppose that it is zero as we work with a small angle fixed to 2.5 deg)
        savedPeakFitParams[1],  ## peak position of HD
        savedPeakFitParams[0],  ## peak intensity of HD (maximum intensity)
        0.5
        * (
            savedPeakFitParams[2] + savedPeakFitParams[3]
        ),  ## FWHM of HD (mean of the FWHM at left and right as I am using an assymetric function)
        savedPeakFitParams[
            4
        ],  ## shape factor (contribution of the lorentzian function)
        savedPeakFitParams[5]
        if len(savedPeakFitParams) >= 6
        else 0,  ## Rw factor of HD (goodness of the fit)
    ]
