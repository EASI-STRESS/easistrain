from typing import Dict, Optional, Sequence, Tuple, Union
import h5py
import numpy


nxchar = h5py.special_dtype(vlen=str)


def as_nxchar(s: Union[str, Sequence[str]]) -> numpy.ndarray:
    return numpy.array(s, dtype=nxchar)


def create_info_group(
    root: h5py.Group,
    fileRead: str,
    fileSave: str,
    sample: Optional[str],
    dataset: Union[str, int, None],
    nameHorizontalDetector: str,
    nameVerticalDetector: str,
    nbPeaksInBoxes: Sequence[int],
):
    infoGroup = root.create_group("infos")  ## infos group creation
    infoGroup.create_dataset(
        "fileRead", dtype=h5py.string_dtype(encoding="utf-8"), data=fileRead
    )  ## save path of raw data file in infos group
    infoGroup.create_dataset(
        "fileSave", dtype=h5py.string_dtype(encoding="utf-8"), data=fileSave
    )  ## save path of the file in which results will be saved in info group
    infoGroup.create_dataset(
        "sample",
        dtype=h5py.string_dtype(encoding="utf-8"),
        data=sample if sample is not None else "No sample",
    )  ## save the name of the sample in infos group
    infoGroup.create_dataset(
        "dataset",
        dtype=h5py.string_dtype(encoding="utf-8"),
        data=dataset if dataset is not None else "No dataset",
    )  ## save the name of dataset in infos group
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
        "nbPeaksInBoxes", dtype="int", data=nbPeaksInBoxes
    )  ## save of the number of peaks per box/window in infos group
    infoGroup.create_dataset(
        "fittingFunction",
        dtype=h5py.string_dtype(encoding="utf-8"),
        data="asymmetric Pseudo-Voigt",
    )  ## save of the type of function used in the fitting of the peaks

    return infoGroup


def create_calib_info_group(
    root: h5py.Group,
    fileRead: str,
    fileSave: str,
    sample: Optional[str],
    dataset: Union[str, int, None],
    nameHorizontalDetector: str,
    nameVerticalDetector: str,
    nbPeaksInBoxes: Sequence[int],
    scanNumberHorizontalDetector: Union[str, int],
    scanNumberVerticalDetector: Union[str, int],
    rangeFit: Sequence[int],
    peakEnergies: Sequence[float],
):
    infoGroup = create_info_group(
        root,
        fileRead,
        fileSave,
        sample,
        dataset,
        nameHorizontalDetector,
        nameVerticalDetector,
        nbPeaksInBoxes,
    )
    infoGroup.create_dataset(
        "scanNumberHorizontalDetector",
        dtype=int,
        data=int(scanNumberHorizontalDetector),
    )  ## save of the number of the scan containing the calibration pattern of the horizontal detector in infos group
    infoGroup.create_dataset(
        "scanNumberVerticalDetector",
        dtype=int,
        data=int(scanNumberVerticalDetector),
    )  ## save of the number of the scan containing the calibration pattern of the vertical detector in info group
    infoGroup.create_dataset(
        "rangeFit", dtype=int, data=rangeFit
    )  ## save of the range of the fit of each box/window in infos group
    infoGroup.create_dataset(
        "peakEnergies",
        dtype=float,
        data=peakEnergies,
    )  ## save of the path of the calibrant File in infos group

    return infoGroup


def create_angle_calib_info_group(
    root: h5py.Group,
    fileRead: str,
    fileSave: str,
    sample: Optional[str],
    dataset: Union[str, int, None],
    nameHorizontalDetector: str,
    nameVerticalDetector: str,
    nbPeaksInBoxes: Sequence[int],
    rangeFitHD: Sequence[int],
    rangeFitVD: Sequence[int],
):
    infoGroup = create_info_group(
        root,
        fileRead,
        fileSave,
        sample,
        dataset,
        nameHorizontalDetector,
        nameVerticalDetector,
        nbPeaksInBoxes,
    )
    infoGroup.create_dataset(
        "rangeFitHD", dtype="int", data=rangeFitHD
    )  ## save of the range of the fit of each box/window of the horizontal detector in infos group
    infoGroup.create_dataset(
        "rangeFitVD", dtype="int", data=rangeFitVD
    )  ## save of the range of the fit of each box/window of the vertical detector in infos group

    return infoGroup


def create_fit_info_group(
    root: h5py.Group,
    fileRead: str,
    fileSave: str,
    sample: Optional[str],
    dataset: Union[str, int, None],
    scanNumber: Union[str, int],
    nameHorizontalDetector: str,
    nameVerticalDetector: str,
    nbPeaksInBoxes: Sequence[int],
    rangeFitHD: Sequence[int],
    rangeFitVD: Sequence[int],
    positioners: Sequence[str],
):
    infoGroup = create_info_group(
        root,
        fileRead,
        fileSave,
        sample,
        dataset,
        nameHorizontalDetector,
        nameVerticalDetector,
        nbPeaksInBoxes,
    )
    infoGroup.create_dataset(
        "scanNumber", dtype="int", data=scanNumber
    )  ## save of the number of the scan in infos group
    infoGroup.create_dataset(
        "rangeFitHD", dtype="int", data=rangeFitHD
    )  ## save of the range of the fit of each box/window of the horizontal detector in infos group
    infoGroup.create_dataset(
        "rangeFitVD", dtype="int", data=rangeFitVD
    )  ## save of the range of the fit of each box/window of the vertical detector in infos group
    infoGroup.create_dataset(
        "positioners", dtype=h5py.string_dtype(encoding="utf-8"), data=str(positioners)
    )  ## save of the range of the fit of each box/window of the vertical detector in infos group

    return infoGroup


def peak_dataset_data(
    positionAngles: numpy.ndarray,
    savedPeakFitParams: numpy.ndarray,
    delta_angle: float,
    nacq,
):
    return [
        *positionAngles[nacq, 0:6],
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


def save_fit_data(
    fitLine: h5py.Group,
    detectorName: str,
    channels: numpy.ndarray,
    raw_data: numpy.ndarray,
    background: numpy.ndarray,
    fitted_data: numpy.ndarray,
):
    detectorGroup = fitLine.create_group(detectorName)
    detectorGroup.create_dataset(
        "channels",
        dtype="float64",
        data=channels,
    )
    detectorGroup.create_dataset(
        "raw_data",
        dtype="float64",
        data=raw_data,
    )
    detectorGroup.create_dataset(
        "background",
        dtype="float64",
        data=background,
    )
    detectorGroup.create_dataset(
        "data - background",
        dtype="float64",
        data=raw_data - background,
    )
    detectorGroup.create_dataset(
        "fitted_data",
        dtype="float64",
        data=fitted_data,
    )
    detectorGroup.create_dataset(
        "residual",
        dtype="float64",
        data=numpy.absolute(fitted_data - raw_data),
    )

    # NeXus
    detectorGroup.attrs["NX_class"] = as_nxchar("NXdata")
    detectorGroup.attrs["auxiliary_signals"] = as_nxchar(["raw_data"])
    detectorGroup.attrs["signal"] = as_nxchar("fitted_data")
    detectorGroup.attrs["axes"] = as_nxchar("channels")


def scan_group_path(
    sample: Union[str, None],
    dataset: Union[str, int, None],
    scanNumber: Union[str, int],
):
    path = ""

    if sample:
        path += f"{sample}_"

    if dataset:
        path += f"{dataset}_"

    return f"{path}{scanNumber}.1"


def read_detector_pattern(
    input_filename: str,
    sample: Union[str, None],
    dataset: Union[str, int, None],
    scanNumber: Union[str, int],
    detector_name: str,
    detectorSliceIndex: Union[int, Tuple[()]] = tuple(),
):
    with h5py.File(input_filename, "r") as h5Read:  ## Read the h5 file of raw data
        meas_group = h5Read[
            f"{scan_group_path(sample, dataset, scanNumber)}/measurement"
        ]
        if not isinstance(meas_group, h5py.Group) or detector_name not in meas_group:
            raise TypeError("No pattern was saved in this scan")
            return
        detector_dset = meas_group[detector_name]
        assert isinstance(detector_dset, h5py.Dataset)
        return detector_dset[detectorSliceIndex]


def read_positioner_data(
    input_filename: str,
    sample: Union[str, None],
    dataset: Union[str, int, None],
    scanNumber: Union[str, int],
    positioner_name: str,
    detectorSliceIndex: Union[int, Tuple[()]] = tuple(),
):
    with h5py.File(input_filename, "r") as h5Read:
        input_positioners = h5Read[
            f"{scan_group_path(sample, dataset, scanNumber)}/instrument/positioners"
        ]
        assert isinstance(input_positioners, h5py.Group)
        if positioner_name not in input_positioners:
            raise KeyError(
                f"{positioner_name} was not found in input positioners! Possible values: {input_positioners.keys()}"
            )
        pos_dataset = input_positioners[positioner_name]
        assert isinstance(pos_dataset, h5py.Dataset)

        # Scalar
        if not pos_dataset.shape:
            return pos_dataset[()]

        return pos_dataset[detectorSliceIndex]


def save_fit_params(
    parent: h5py.Group,
    fitParams: Dict[str, numpy.ndarray],
    uncertaintyFitParams: Dict[str, numpy.ndarray],
):
    savedFitParamsHD = numpy.reshape(
        fitParams["horizontal"], (int(numpy.size(fitParams["horizontal"]) / 6), 6)
    )
    parent.create_dataset(
        "fitParamsHD",
        dtype="float64",
        data=savedFitParamsHD,
    )  ## save parameters of the fit of HD
    savedUncertaintyFitParamsHD = numpy.reshape(
        uncertaintyFitParams["horizontal"],
        (int(numpy.size(uncertaintyFitParams["horizontal"]) / 5), 5),
    )
    parent.create_dataset(
        "uncertaintyFitParamsHD",
        dtype="float64",
        data=savedUncertaintyFitParamsHD,
    )  ## save uncertainty on the parameters of the fit of HD

    savedFitParamsVD = numpy.reshape(
        fitParams["vertical"], (int(numpy.size(fitParams["vertical"]) / 6), 6)
    )
    parent.create_dataset(
        "fitParamsVD",
        dtype="float64",
        data=savedFitParamsVD,
    )  ## save parameters of the fit of VD
    savedUncertaintyFitParamsVD = numpy.reshape(
        uncertaintyFitParams["vertical"],
        (int(numpy.size(uncertaintyFitParams["vertical"]) / 5), 5),
    )
    parent.create_dataset(
        "uncertaintyFitParamsVD",
        dtype="float64",
        data=savedUncertaintyFitParamsVD,
    )  ## save uncertainty on the parameters of the fit of VD

    return (
        savedFitParamsHD,
        savedFitParamsVD,
        savedUncertaintyFitParamsHD,
        savedUncertaintyFitParamsVD,
    )
