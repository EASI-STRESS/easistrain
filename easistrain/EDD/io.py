from typing import Sequence, Union
import h5py
import numpy as np

nxchar = h5py.special_dtype(vlen=str)


def as_nxchar(s: Union[str, Sequence[str]]) -> np.ndarray:
    return np.array(s, dtype=nxchar)


def create_info_group(
    root: h5py.Group,
    fileRead: str,
    fileSave: str,
    sample: str,
    dataset: str,
    scanNumber: int,
    nameHorizontalDetector: str,
    nameVerticalDetector: str,
    numberOfBoxes: int,
    nbPeaksInBoxes: Sequence[int],
    rangeFitHD: Sequence[int],
    rangeFitVD: Sequence[int],
    positioners: Sequence[str],
):
    infoGroup = root.create_group("infos")  ## infos group creation
    infoGroup.create_dataset(
        "fileRead", dtype=h5py.string_dtype(encoding="utf-8"), data=fileRead
    )  ## save path of raw data file in infos group
    infoGroup.create_dataset(
        "fileSave", dtype=h5py.string_dtype(encoding="utf-8"), data=fileSave
    )  ## save path of the file in which results will be saved in info group
    infoGroup.create_dataset(
        "sample", dtype=h5py.string_dtype(encoding="utf-8"), data=sample
    )  ## save the name of the sample in infos group
    infoGroup.create_dataset(
        "dataset", dtype=h5py.string_dtype(encoding="utf-8"), data=dataset
    )  ## save the name of dataset in infos group
    infoGroup.create_dataset(
        "scanNumber", dtype="int", data=scanNumber
    )  ## save of the number of the scan in infos group
    infoGroup.create_dataset(
        "nameHorizontalDetector",
        dtype=h5py.string_dtype(encoding="utf-8"),
        data=nameHorizontalDetector,
    )  ## save of the name of the horizontal detector in infos group
    infoGroup.create_dataset(
        "nameVerticalDetector",
        dtype=h5py.string_dtype(encoding="utf-8"),
        data=nameVerticalDetector,
    )  ## save of the name of the vertical detector in infos group
    infoGroup.create_dataset(
        "numberOfBoxes", dtype="int", data=numberOfBoxes
    )  ## save of the number of the boxes/widows extracted from the raw data in infos group
    infoGroup.create_dataset(
        "nbPeaksInBoxes", dtype="int", data=nbPeaksInBoxes
    )  ## save of the number of peaks per box/window in infos group
    infoGroup.create_dataset(
        "rangeFitHD", dtype="int", data=rangeFitHD
    )  ## save of the range of the fit of each box/window of the horizontal detector in infos group
    infoGroup.create_dataset(
        "rangeFitVD", dtype="int", data=rangeFitVD
    )  ## save of the range of the fit of each box/window of the vertical detector in infos group
    infoGroup.create_dataset(
        "positioners", dtype=h5py.string_dtype(encoding="utf-8"), data=str(positioners)
    )  ## save of the range of the fit of each box/window of the vertical detector in infos group


def peak_dataset_data(
    positionAngles: np.ndarray, savedPeakFitParams: np.ndarray, delta_angle: float
):
    return [
        *positionAngles[0, 0:6],
        delta_angle,  ## delta angle of the horizontal detector (debye scherer angle)
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
