import h5py
import numpy
from typing import Sequence, Union
from easistrain.EDD.detector_fit import fit_all_peaks_and_save_results
from easistrain.EDD.io import create_calib_info_group, read_detector_pattern
from easistrain.EDD.utils import run_from_cli
from easistrain.calibrants import calibrant_filename


def calibEdd(
    fileRead: str,
    fileSave: str,
    sample: str,
    dataset: Union[str, int],
    scanNumberHorizontalDetector: Union[str, int],
    scanNumberVerticalDetector: Union[str, int],
    nameHorizontalDetector: str,
    nameVerticalDetector: str,
    nbPeaksInBoxes: Sequence[int],
    rangeFit: Sequence[int],
    sourceCalibrantFile: str,
):
    """Energy calibration of the channels of an energy-dispersive detector (EDD)"""

    patternHorizontalDetector = read_detector_pattern(
        fileRead, sample, dataset, scanNumberHorizontalDetector, nameHorizontalDetector
    )[0]
    patternVerticalDetector = read_detector_pattern(
        fileRead, sample, dataset, scanNumberVerticalDetector, nameVerticalDetector
    )[0]

    with h5py.File(fileSave, "a") as h5Save:
        if "detectorCalibration" not in h5Save.keys():
            calibrationLevel1 = h5Save.create_group(
                "detectorCalibration"
            )  ## calibration group
        else:
            calibrationLevel1 = h5Save["detectorCalibration"]
            assert isinstance(calibrationLevel1, h5py.Group)
        rawDataLevel1_1 = calibrationLevel1.create_group(
            f"rawData_{dataset}_{scanNumberHorizontalDetector}_{scanNumberVerticalDetector}"
        )  ## rawData subgroup in calibration group
        fitLevel1_2 = calibrationLevel1.create_group(
            f"fit_{dataset}_{scanNumberHorizontalDetector}_{scanNumberVerticalDetector}"
        )  ## fit subgroup in calibration group
        fitLevel1_2.create_group(
            "curveCalibration"
        )  ## curve calibration group for the two detector
        fitLevel1_2.create_group(
            "calibCoeffs"
        )  ## calibration coefficients group for the two detector

        create_calib_info_group(
            fitLevel1_2,
            fileRead,
            fileSave,
            sample,
            dataset,
            nameHorizontalDetector,
            nameVerticalDetector,
            nbPeaksInBoxes,
            scanNumberHorizontalDetector,
            scanNumberVerticalDetector,
            rangeFit,
            sourceCalibrantFile,
        )

        curveCalibrationHD = numpy.zeros((numpy.sum(nbPeaksInBoxes), 2), float)
        curveCalibrationVD = numpy.zeros((numpy.sum(nbPeaksInBoxes), 2), float)

        fit_all_peaks_and_save_results(
            nbPeaksInBoxes,
            rangeFit={"horizontal": rangeFit, "vertical": rangeFit},
            patterns={
                "horizontal": patternHorizontalDetector,
                "vertical": patternVerticalDetector,
            },
            scanNumbers={
                "horizontal": scanNumberHorizontalDetector,
                "vertical": scanNumberVerticalDetector,
            },
            saving_dest=fitLevel1_2,
            group_format=lambda i: f"fitLine_{i}",
        )

        rawDataLevel1_1.create_dataset(
            "horizontalDetector", dtype="float64", data=patternHorizontalDetector
        )  ## save raw data of the horizontal detector
        rawDataLevel1_1.create_dataset(
            "verticalDetector", dtype="float64", data=patternVerticalDetector
        )  ## save raw data of the vertical detector

        calibrantSource = numpy.loadtxt(
            calibrant_filename(sourceCalibrantFile)
        )  ## open source calibration text file
        curveCalibrationHD[:, 0] = fitLevel1_2["fitParams/fitParamsHD"][:, 1]
        curveCalibrationHD[:, 1] = calibrantSource[: numpy.sum(nbPeaksInBoxes)]
        ucurveCalibrationHD = fitLevel1_2["fitParams/uncertaintyFitParamsHD"][:, 1]
        fitLevel1_2["curveCalibration"].create_dataset(
            "curveCalibrationHD", dtype="float64", data=curveCalibrationHD
        )  ## curve energy VS channels for horizontal detector
        curveCalibrationVD[:, 0] = fitLevel1_2["fitParams/fitParamsVD"][:, 1]
        curveCalibrationVD[:, 1] = calibrantSource[: numpy.sum(nbPeaksInBoxes)]
        ucurveCalibrationVD = fitLevel1_2["fitParams/uncertaintyFitParamsVD"][:, 1]
        fitLevel1_2["curveCalibration"].create_dataset(
            "curveCalibrationVD", dtype="float64", data=curveCalibrationVD
        )  ## curve energy VS channels for vertical detector
        calibCoeffsHD, covCalibCoeffsHD = numpy.polyfit(
            x=curveCalibrationHD[:, 0],
            y=curveCalibrationHD[:, 1],
            deg=2,
            full=False,
            cov="unscaled",
            w=1 / ucurveCalibrationHD,
        )  ## calibration coefficients of the horizontal detector
        calibCoeffsVD, covCalibCoeffsVD = numpy.polyfit(
            x=curveCalibrationVD[:, 0],
            y=curveCalibrationVD[:, 1],
            deg=2,
            full=False,
            cov="unscaled",
            w=1 / ucurveCalibrationVD,
        )  ## calibration coefficients of the vertical detector
        fitLevel1_2["curveCalibration"].create_dataset(
            "fitCurveCalibrationHD",
            dtype="float64",
            data=numpy.transpose(
                (
                    curveCalibrationHD[:, 0],
                    numpy.poly1d(calibCoeffsHD)(curveCalibrationHD[:, 0]),
                )
            ),
        )  ## fitted curve energy VS channels for horizontal detector
        fitLevel1_2["curveCalibration"].create_dataset(
            "fitCurveCalibrationVD",
            dtype="float64",
            data=numpy.transpose(
                (
                    curveCalibrationVD[:, 0],
                    numpy.poly1d(calibCoeffsVD)(curveCalibrationVD[:, 0]),
                )
            ),
        )  ## fitted curve energy VS channels for vertical detector
        fitLevel1_2["curveCalibration"].create_dataset(
            "errorCurveCalibrationHD",
            dtype="float64",
            data=numpy.transpose(
                (
                    curveCalibrationHD[:, 0],
                    numpy.abs(
                        numpy.poly1d(calibCoeffsHD)(curveCalibrationHD[:, 0])
                        - curveCalibrationHD[:, 1]
                    ),
                )
            ),
        )  ## error between fitted and raw curve energy VS channels for horizontal detector
        fitLevel1_2["curveCalibration"].create_dataset(
            "errorCurveCalibrationVD",
            dtype="float64",
            data=numpy.transpose(
                (
                    curveCalibrationVD[:, 0],
                    numpy.abs(
                        numpy.poly1d(calibCoeffsVD)(curveCalibrationVD[:, 0])
                        - curveCalibrationVD[:, 1]
                    ),
                )
            ),
        )  ## error between fitted and raw curve energy VS channels for vertical detector
        # print(f'uncertauntyCalibCoeffsHD = {numpy.sqrt(numpy.diag(covCalibCoeffsHD))}')
        # print(f'uncertauntyCalibCoeffsVD = {numpy.sqrt(numpy.diag(covCalibCoeffsVD))}')
        fitLevel1_2["calibCoeffs"].create_dataset(
            "calibCoeffsHD", dtype="float64", data=calibCoeffsHD
        )  ## save calibration coefficients of the horizontal detector
        fitLevel1_2["calibCoeffs"].create_dataset(
            "calibCoeffsVD", dtype="float64", data=calibCoeffsVD
        )  ## save calibration coefficients of the Vertical detector
        fitLevel1_2["calibCoeffs"].create_dataset(
            "uncertaintyCalibCoeffsHD",
            dtype="float64",
            data=numpy.sqrt(numpy.diag(covCalibCoeffsHD)),
        )  ## save uncertainty on calibration coefficients of the horizontal detector
        fitLevel1_2["calibCoeffs"].create_dataset(
            "uncertaintyCalibCoeffsVD",
            dtype="float64",
            data=numpy.sqrt(numpy.diag(covCalibCoeffsVD)),
        )  ## save uncertainty on calibration coefficients of the vertical detector


if __name__ == "__main__":
    run_from_cli(calibEdd)
