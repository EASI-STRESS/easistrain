import os
from typing import Optional, Sequence, Union

import h5py
import numpy
import scipy.optimize
import scipy.constants

from easistrain.EDD.constants import pCstInkeVS, speedLightInAPerS
from easistrain.EDD.detector_fit import fit_all_peaks_and_save_results
from easistrain.EDD.io import create_angle_calib_info_group, read_detector_pattern
from easistrain.EDD.utils import linefunc, run_from_cli, uChEConversion
from easistrain.calibrants import calibrant_filename


def angleCalibrationEDD(
    fileRead: str,
    fileSave: str,
    sample: Optional[str],
    dataset: Optional[str],
    scanNumber: Union[str, int],
    nameHorizontalDetector: str,
    nameVerticalDetector: str,
    nbPeaksInBoxes: Sequence[int],
    rangeFitHD: Sequence[int],
    rangeFitVD: Sequence[int],
    pathFileDetectorCalibration: str,
    scanDetectorCalibration: str,
    sampleCalibrantFile: str,
):
    """Main function."""
    nbPeaks = sum(nbPeaksInBoxes)

    patternHorizontalDetector = read_detector_pattern(
        fileRead, sample, dataset, scanNumber, nameHorizontalDetector
    )
    patternVerticalDetector = read_detector_pattern(
        fileRead, sample, dataset, scanNumber, nameVerticalDetector
    )
    if patternHorizontalDetector.ndim == 2:
        patternHorizontalDetector = patternHorizontalDetector[0]
    if patternVerticalDetector.ndim == 2:
        patternVerticalDetector = patternVerticalDetector[0]

    if os.path.dirname(fileSave):
        os.makedirs(os.path.dirname(fileSave), exist_ok=True)
    with h5py.File(fileSave, "a") as h5Save:  ## create/append h5 file to save in
        if "angleCalibration" not in h5Save.keys():
            angleCalibrationLevel1 = h5Save.create_group(
                "angleCalibration"
            )  ## angleCalibration group
        else:
            angleCalibrationLevel1 = h5Save["angleCalibration"]
        rawDataLevel1_1 = angleCalibrationLevel1.create_group(
            "_".join(
                [str(v) for v in ["rawData", dataset, scanNumber] if v is not None]
            )
        )  ## rawData subgroup in calibration group
        fitLevel1_2 = angleCalibrationLevel1.create_group(
            "_".join([str(v) for v in ["fit", dataset, scanNumber] if v is not None])
        )  ## fit subgroup in calibration group
        fitLevel1_2.create_group(
            "curveAngleCalibration"
        )  ## curve E VS channels group for the two detector
        fitLevel1_2.create_group(
            "calibratedAngle"
        )  ## diffraction angle of the two detectors group for the two detector

        infoGroup = create_angle_calib_info_group(
            fitLevel1_2,
            fileRead,
            fileSave,
            sample,
            dataset,
            nameHorizontalDetector,
            nameVerticalDetector,
            nbPeaksInBoxes,
            rangeFitHD,
            rangeFitVD,
        )

        curveAngleCalibrationHD = numpy.zeros((numpy.sum(nbPeaksInBoxes), 2), float)
        curveAngleCalibrationVD = numpy.zeros((numpy.sum(nbPeaksInBoxes), 2), float)

        (savedFitParamsHD, savedFitParamsVD, _, __) = fit_all_peaks_and_save_results(
            nbPeaksInBoxes,
            rangeFit={"horizontal": rangeFitHD, "vertical": rangeFitVD},
            patterns={
                "horizontal": patternHorizontalDetector,
                "vertical": patternVerticalDetector,
            },
            scanNumbers={
                "horizontal": scanNumber,
                "vertical": scanNumber,
            },
            saving_dest=fitLevel1_2,
            group_format=lambda i: f"fitLine_{str(i).zfill(4)}",
        )

        rawDataLevel1_1.create_dataset(
            "horizontalDetector", dtype="float64", data=patternHorizontalDetector
        )  ## save raw data of the horizontal detector
        rawDataLevel1_1.create_dataset(
            "verticalDetector", dtype="float64", data=patternVerticalDetector
        )  ## save raw data of the vertical detector

        with h5py.File(
            pathFileDetectorCalibration, "r"
        ) as h5pFDC:  ## Read the h5 file of the energy calibration of the two detectors
            calibCoeffsHD = h5pFDC[
                f"detectorCalibration/{scanDetectorCalibration}/calibCoeffs/calibCoeffsHD"
            ][
                ()
            ]  ## import the energy calibration coefficients of the horizontal detector
            calibCoeffsVD = h5pFDC[
                f"detectorCalibration/{scanDetectorCalibration}/calibCoeffs/calibCoeffsVD"
            ][
                ()
            ]  ## import the energy calibration coefficients of the vertical detector
            uncertaintyCalibCoeffsHD = h5pFDC[
                f"detectorCalibration/{scanDetectorCalibration}/calibCoeffs/uncertaintyCalibCoeffsHD"
            ][
                ()
            ]  ## import the uncertainty of the energy calibration coefficients of the horizontal detector
            uncertaintyCalibCoeffsVD = h5pFDC[
                f"detectorCalibration/{scanDetectorCalibration}/calibCoeffs/uncertaintyCalibCoeffsVD"
            ][
                ()
            ]  ## import the uncertainty of the energy calibration coefficients of the vertical detector
        calibrantSample = numpy.loadtxt(
            calibrant_filename(sampleCalibrantFile)
        )  ## open source calibration text file
        if len(calibrantSample) < nbPeaks:
            raise ValueError(
                f"d-spacing file {sampleCalibrantFile} should have at leasts {nbPeaks} peaks"
            )

        conversionChannelEnergyHD = numpy.polyval(
            calibCoeffsHD, savedFitParamsHD[:, 1]
        )  ## conversion of the channel to energy for the horizontal detector
        conversionChannelEnergyVD = numpy.polyval(
            calibCoeffsVD, savedFitParamsVD[:, 1]
        )  ## conversion of the channel to energy for the vertical detector
        curveAngleCalibrationHD[:, 0] = 1 / calibrantSample[:nbPeaks]
        curveAngleCalibrationHD[:, 1] = conversionChannelEnergyHD
        fitLevel1_2["curveAngleCalibration"].create_dataset(
            "curveAngleCalibrationHD", dtype="float64", data=curveAngleCalibrationHD
        )  ## save curve energy VS 1/d for horizontal detector (d = hkl interriticular distance of the calibrant sample)
        curveAngleCalibrationVD[:, 0] = 1 / calibrantSample[:nbPeaks]
        curveAngleCalibrationVD[:, 1] = conversionChannelEnergyVD
        fitLevel1_2["curveAngleCalibration"].create_dataset(
            "curveAngleCalibrationVD", dtype="float64", data=curveAngleCalibrationVD
        )  ## save curve energy VS 1/d for vertical detector (d = hkl interriticular distance of the calibrant sample)

        calibratedAngleHD, covCalibratedAngleHD = scipy.optimize.curve_fit(
            f=linefunc,
            xdata=curveAngleCalibrationHD[:, 0],
            ydata=curveAngleCalibrationHD[:, 1],
            p0=numpy.polyfit(
                x=curveAngleCalibrationHD[:, 0], y=curveAngleCalibrationHD[:, 1], deg=1
            )[0],
            sigma=uChEConversion(
                calibCoeffsHD[0],
                calibCoeffsHD[1],
                calibCoeffsHD[2],
                fitLevel1_2["fitParams/fitParamsHD"][:, 1],
                uncertaintyCalibCoeffsHD[0],
                uncertaintyCalibCoeffsHD[1],
                uncertaintyCalibCoeffsHD[2],
                fitLevel1_2["fitParams/uncertaintyFitParamsHD"][:, 1],
            ),
        )  ## calculation of 12.398/2*sin(theta) of the diffraction angle of the horizontal detector
        calibratedAngleVD, covCalibratedAngleVD = scipy.optimize.curve_fit(
            f=linefunc,
            xdata=curveAngleCalibrationVD[:, 0],
            ydata=curveAngleCalibrationVD[:, 1],
            p0=numpy.polyfit(
                x=curveAngleCalibrationVD[:, 0], y=curveAngleCalibrationVD[:, 1], deg=1
            )[0],
            sigma=uChEConversion(
                calibCoeffsVD[0],
                calibCoeffsVD[1],
                calibCoeffsVD[2],
                fitLevel1_2["fitParams/fitParamsHD"][:, 1],
                uncertaintyCalibCoeffsVD[0],
                uncertaintyCalibCoeffsVD[1],
                uncertaintyCalibCoeffsVD[2],
                fitLevel1_2["fitParams/uncertaintyFitParamsVD"][:, 1],
            ),
        )  ## calculation of 12.398/2*sin(theta) of the diffraction angle of the vertical detector
        # print(calibratedAngleHD)
        # print(calibratedAngleVD)
        fitLevel1_2["calibratedAngle"].create_dataset(
            "calibratedAngleHD",
            dtype="float64",
            data=numpy.rad2deg(
                2
                * numpy.arcsin(
                    (pCstInkeVS * speedLightInAPerS) / (2 * calibratedAngleHD)
                )
            ),
        )  ## save the calibrated diffraction angle in degree of the horizontal detector
        fitLevel1_2["calibratedAngle"].create_dataset(
            "calibratedAngleVD",
            dtype="float64",
            data=numpy.rad2deg(
                2
                * numpy.arcsin(
                    (pCstInkeVS * speedLightInAPerS) / (2 * calibratedAngleVD)
                )
            ),
        )  ## save the calibrated diffraction angle in degree of the vertical detector
        fitLevel1_2["calibratedAngle"].create_dataset(
            "uncertaintyCalibratedAngleHD",
            dtype="float64",
            data=numpy.sqrt(numpy.diag(covCalibratedAngleHD)),
        )  ## save the uncertainty of the calibrated diffraction angle in degree of the horizontal detector
        fitLevel1_2["calibratedAngle"].create_dataset(
            "uncertaintyCalibratedAngleVD",
            dtype="float64",
            data=numpy.sqrt(numpy.diag(covCalibratedAngleVD)),
        )  ## save the uncertainty of the calibrated diffraction angle in degree of the vertical detector
        fitLevel1_2["curveAngleCalibration"].create_dataset(
            "fitCurveAngleCalibrationHD",
            dtype="float64",
            data=numpy.transpose(
                (
                    curveAngleCalibrationHD[:, 0],
                    calibratedAngleHD * curveAngleCalibrationHD[:, 0],
                )
            ),
        )  ## save curve energy VS 1/d for horizontal detector calculated using the fitted value 12.398/2*sin(theta)
        fitLevel1_2["curveAngleCalibration"].create_dataset(
            "fitCurveAngleCalibrationVD",
            dtype="float64",
            data=numpy.transpose(
                (
                    curveAngleCalibrationVD[:, 0],
                    calibratedAngleVD * curveAngleCalibrationVD[:, 0],
                )
            ),
        )  ## save curve energy VS 1/d for vertical detector calculated using the fitted value 12.398/2*sin(theta)
        fitLevel1_2["curveAngleCalibration"].create_dataset(
            "errorCurveAngleCalibrationHD",
            dtype="float64",
            data=numpy.transpose(
                (
                    curveAngleCalibrationHD[:, 0],
                    numpy.abs(
                        curveAngleCalibrationHD[:, 1]
                        - (calibratedAngleHD * curveAngleCalibrationHD[:, 0])
                    ),
                )
            ),
        )  ## error between the fitted (calculated using the fitted value 12.398/2*sin(theta)) and experimental curve of energy VS 1/d for horizontal detector
        fitLevel1_2["curveAngleCalibration"].create_dataset(
            "errorCurveAngleCalibrationVD",
            dtype="float64",
            data=numpy.transpose(
                (
                    curveAngleCalibrationVD[:, 0],
                    numpy.abs(
                        curveAngleCalibrationVD[:, 1]
                        - (calibratedAngleVD * curveAngleCalibrationVD[:, 0])
                    ),
                )
            ),
        )  ## error between the fitted (calculated using the fitted value 12.398/2*sin(theta)) and experimental curve of energy VS 1/d for vertical detector

        infoGroup.create_dataset(
            "pathDetectorCalibrationParams",
            dtype=h5py.string_dtype(encoding="utf-8"),
            data=f"{pathFileDetectorCalibration}/detectorCalibration/{scanDetectorCalibration}/calibCoeffs",
        )  ## save of the path of the file containing the energy calibration coefficient for the two detectors used for the conversion of channels ====> energy in the info group


if __name__ == "__main__":
    run_from_cli(angleCalibrationEDD)
