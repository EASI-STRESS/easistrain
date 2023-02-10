import h5py
import numpy as np
import silx.math.fit
import silx.math.fit.peaks
import scipy.optimize
from typing import Sequence, Union
from easistrain.EDD.io import create_calib_info_group

from easistrain.EDD.utils import (
    calcBackground,
    guessParameters,
    run_from_cli,
    splitPseudoVoigt,
)


def calibEdd(
    fileRead: str,
    fileSave: str,
    sample: str,
    dataset: Union[str, int],
    scanNumberHorizontalDetector: Union[str, int],
    scanNumberVerticalDetector: Union[str, int],
    nameHorizontalDetector: str,
    nameVerticalDetector: str,
    numberOfBoxes: int,
    nbPeaksInBoxes: Sequence[int],
    rangeFit: Sequence[int],
    sourceCalibrantFile: str,
):
    """Main function."""

    with h5py.File(fileRead, "r") as h5Read:  ## Read the h5 file of raw data
        patternHorizontalDetector = h5Read[
            f"{sample}_{dataset}_{scanNumberHorizontalDetector}.1/measurement/{nameHorizontalDetector}"
        ][0]
        patternVerticalDetector = h5Read[
            f"{sample}_{dataset}_{scanNumberVerticalDetector}.1/measurement/{nameVerticalDetector}"
        ][0]

    h5Save = h5py.File(fileSave, "a")  ## create h5 file to save in
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
    fitLevel1_2.create_group("fitParams")  ## fit results group for the two detector
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
        numberOfBoxes,
        nbPeaksInBoxes,
        scanNumberHorizontalDetector,
        scanNumberVerticalDetector,
        rangeFit,
        sourceCalibrantFile,
    )

    fitParamsHD = np.array(())
    fitParamsVD = np.array(())
    uncertaintyFitParamsHD = np.array(())
    uncertaintyFitParamsVD = np.array(())
    curveCalibrationHD = np.zeros((np.sum(nbPeaksInBoxes), 2), float)
    curveCalibrationVD = np.zeros((np.sum(nbPeaksInBoxes), 2), float)
    for i in range(numberOfBoxes):
        peakHorizontalDetector = np.transpose(
            (
                np.arange(rangeFit[2 * i], rangeFit[(2 * i) + 1]),
                patternHorizontalDetector[rangeFit[2 * i] : rangeFit[(2 * i) + 1]],
            )
        )  ## peak of the horizontal detector
        peakVerticalDetector = np.transpose(
            (
                np.arange(rangeFit[2 * i], rangeFit[(2 * i) + 1]),
                patternVerticalDetector[rangeFit[2 * i] : rangeFit[(2 * i) + 1]],
            )
        )  ## peak of the vertical detector
        backgroundHorizontalDetector = silx.math.fit.strip(
            data=peakHorizontalDetector[:, 1],
            w=5,
            niterations=5000,
            factor=1,
            anchors=None,
        )  ## stripped background of the horizontal detector (obtained by stripping the yData)
        backgroundVerticalDetector = silx.math.fit.strip(
            data=peakVerticalDetector[:, 1],
            w=5,
            niterations=5000,
            factor=1,
            anchors=None,
        )  ## stripped background of the vertical detector (obtained by stripping the yData)
        # print(backgroundHorizontalDetector)
        # print(backgroundVerticalDetector)
        fit_line_group = fitLevel1_2.create_group(
            f"fitLine_{i}"
        )  ## create group for each calibration peak
        fit_line_group.create_dataset(
            "rawHorizontalDetector", dtype="float64", data=peakHorizontalDetector
        )  ## create dataset for raw data of each calibration peak
        fit_line_group.create_dataset(
            "rawVerticalDetector", dtype="float64", data=peakVerticalDetector
        )  ## create dataset for raw data of each calibration peak

        peaksGuessHD, peaksIndexHD = guessParameters(
            peakHorizontalDetector[:, 0],
            peakHorizontalDetector[:, 1] - backgroundHorizontalDetector,
            nbPeaksInBoxes[i],
            withBounds=False,
        )  ## guess fit parameters for HD
        # print(peaksIndexHD)
        peaksGuessVD, peaksIndexVD = guessParameters(
            peakVerticalDetector[:, 0],
            peakVerticalDetector[:, 1] - backgroundVerticalDetector,
            nbPeaksInBoxes[i],
            withBounds=False,
        )  ## guess fit parameters for VD
        # print(peaksIndexVD)
        yCalculatedBackgroundHD = calcBackground(
            peakHorizontalDetector[:, 0],
            peakHorizontalDetector[:, 1],
            peaksGuessHD[-1],
            peaksGuessHD[2],
            peaksIndexHD,
        )  ## calculated ybackground of the horizontal detector
        yCalculatedBackgroundVD = calcBackground(
            peakVerticalDetector[:, 0],
            peakVerticalDetector[:, 1],
            peaksGuessVD[-1],
            peaksGuessVD[2],
            peaksIndexVD,
        )  ## calculated ybackground of the vertical detector
        fit_line_group.create_dataset(
            "backgroundHorizontalDetector",
            dtype="float64",
            data=np.transpose((peakHorizontalDetector[:, 0], yCalculatedBackgroundHD)),
        )  ## create dataset for background of each calibration peak for HD
        fit_line_group.create_dataset(
            "backgroundVerticalDetector",
            dtype="float64",
            data=np.transpose((peakVerticalDetector[:, 0], yCalculatedBackgroundVD)),
        )  ## create dataset for background of each calibration peak for VD
        fit_line_group.create_dataset(
            "bgdSubsDataHorizontalDetector",
            dtype="float64",
            data=np.transpose(
                (
                    peakHorizontalDetector[:, 0],
                    peakHorizontalDetector[:, 1] - yCalculatedBackgroundHD,
                )
            ),
        )  ## create dataset for HD raw data after subst of background
        fit_line_group.create_dataset(
            "bgdSubsDataVerticalDetector",
            dtype="float64",
            data=np.transpose(
                (
                    peakVerticalDetector[:, 0],
                    peakVerticalDetector[:, 1] - yCalculatedBackgroundVD,
                )
            ),
        )  ## create dataset for VD raw data after subst of background
        # print(peaksGuessHD)
        # print(peaksGuessVD)
        initialGuessHD = np.zeros(5 * nbPeaksInBoxes[i])
        initialGuessVD = np.zeros(5 * nbPeaksInBoxes[i])
        # print(peaksGuessHD, peaksGuessVD)
        # print(([0,
        #    np.amin(peakHorizontalDetector[:, 0]), 0, 0, 0],
        #   [np.amax(peakHorizontalDetector[:, 1]),
        #  np.amax(peakHorizontalDetector[:, 0]), len(peakHorizontalDetector[:, 0]),
        #  len(peakHorizontalDetector[:, 0]), 1]))
        fit_min_boundsHD = np.zeros(5 * nbPeaksInBoxes[i])
        fit_max_boundsHD = np.zeros(5 * nbPeaksInBoxes[i])
        fit_min_boundsVD = np.zeros(5 * nbPeaksInBoxes[i])
        fit_max_boundsVD = np.zeros(5 * nbPeaksInBoxes[i])
        for n in range(nbPeaksInBoxes[i]):
            initialGuessHD[5 * n] = peaksGuessHD[3 * n]
            initialGuessHD[5 * n + 1] = peaksGuessHD[3 * n + 1]
            initialGuessHD[5 * n + 2] = peaksGuessHD[3 * n + 2]
            initialGuessHD[5 * n + 3] = peaksGuessHD[3 * n + 2]
            initialGuessHD[5 * n + 4] = 0.5
            initialGuessVD[5 * n] = peaksGuessVD[3 * n]
            initialGuessVD[5 * n + 1] = peaksGuessVD[3 * n + 1]
            initialGuessVD[5 * n + 2] = peaksGuessVD[3 * n + 2]
            initialGuessVD[5 * n + 3] = peaksGuessVD[3 * n + 2]
            initialGuessVD[5 * n + 4] = 0.5
            fit_min_boundsHD[5 * n : 5 * n + 5] = [
                0,
                np.amin(peakHorizontalDetector[:, 0]),
                0,
                0,
                0,
            ]
            fit_max_boundsHD[5 * n : 5 * n + 5] = [
                np.inf,
                np.amax(peakHorizontalDetector[:, 0]),
                len(peakHorizontalDetector[:, 0]) / 2,
                len(peakHorizontalDetector[:, 0]) / 2,
                1,
            ]
            fit_min_boundsVD[5 * n : 5 * n + 5] = [
                0,
                np.amin(peakVerticalDetector[:, 0]),
                0,
                0,
                0,
            ]
            fit_max_boundsVD[5 * n : 5 * n + 5] = [
                np.inf,
                np.amax(peakVerticalDetector[:, 0]),
                len(peakVerticalDetector[:, 0]) / 2,
                len(peakVerticalDetector[:, 0]) / 2,
                1,
            ]
        optimal_parametersHD, covarianceHD = scipy.optimize.curve_fit(
            f=splitPseudoVoigt,
            xdata=peakHorizontalDetector[:, 0],
            ydata=peakHorizontalDetector[:, 1] - yCalculatedBackgroundHD,
            p0=initialGuessHD,
            sigma=np.sqrt(0.5 + peakHorizontalDetector[:, 1]),
            bounds=(fit_min_boundsHD, fit_max_boundsHD),
        )  ## fit of the peak of the Horizontal detector
        optimal_parametersVD, covarianceVD = scipy.optimize.curve_fit(
            f=splitPseudoVoigt,
            xdata=peakVerticalDetector[:, 0],
            ydata=peakVerticalDetector[:, 1] - yCalculatedBackgroundVD,
            p0=initialGuessVD,
            sigma=np.sqrt(0.5 + peakVerticalDetector[:, 1]),
            bounds=(fit_min_boundsVD, fit_max_boundsVD),
        )  ## fit of the peak of the Vertical detector
        fit_line_group.create_dataset(
            "fitHorizontalDetector",
            dtype="float64",
            data=np.transpose(
                (
                    peakHorizontalDetector[:, 0],
                    splitPseudoVoigt(peakHorizontalDetector[:, 0], optimal_parametersHD)
                    + yCalculatedBackgroundHD,
                )
            ),
        )  ## fitted data of the horizontal detector
        fit_line_group.create_dataset(
            "fitVerticalDetector",
            dtype="float64",
            data=np.transpose(
                (
                    peakVerticalDetector[:, 0],
                    splitPseudoVoigt(peakVerticalDetector[:, 0], optimal_parametersVD)
                    + yCalculatedBackgroundVD,
                )
            ),
        )  ## fitted data of the vertical detector
        fit_line_group.create_dataset(
            "errorHorizontalDetector",
            dtype="float64",
            data=np.transpose(
                (
                    peakHorizontalDetector[:, 0],
                    np.absolute(
                        splitPseudoVoigt(
                            peakHorizontalDetector[:, 0], optimal_parametersHD
                        )
                        + yCalculatedBackgroundHD
                        - peakHorizontalDetector[:, 1]
                    ),
                )
            ),
        )  ## error of the horizontal detector
        fit_line_group.create_dataset(
            "errorVerticalDetector",
            dtype="float64",
            data=np.transpose(
                (
                    peakVerticalDetector[:, 0],
                    np.absolute(
                        splitPseudoVoigt(
                            peakVerticalDetector[:, 0], optimal_parametersVD
                        )
                        + yCalculatedBackgroundVD
                        - peakVerticalDetector[:, 1]
                    ),
                )
            ),
        )  ## error of the vertical detector
        # print(f'optimal_parametersHD = {optimal_parametersHD}')
        # print(f'uncertauntyHD = {np.sqrt(np.diag(covarianceHD))}')
        # print(f'optimal_parametersVD = {optimal_parametersVD}')
        # print(f'uncertauntyVD = {np.sqrt(np.diag(covarianceVD))}')
        for n in range(nbPeaksInBoxes[i]):
            fitParamsHD = np.append(
                fitParamsHD,
                np.append(
                    optimal_parametersHD[5 * n : 5 * n + 5],
                    100
                    * np.sum(
                        np.absolute(
                            splitPseudoVoigt(
                                peakHorizontalDetector[:, 0], optimal_parametersHD
                            )
                            + backgroundHorizontalDetector
                            - peakHorizontalDetector[:, 1]
                        )
                    )
                    / np.sum(peakHorizontalDetector[:, 1]),
                ),
                axis=0,
            )  ##
            fitParamsVD = np.append(
                fitParamsVD,
                np.append(
                    optimal_parametersVD[5 * n : 5 * n + 5],
                    100
                    * np.sum(
                        np.absolute(
                            splitPseudoVoigt(
                                peakVerticalDetector[:, 0], optimal_parametersVD
                            )
                            + backgroundVerticalDetector
                            - peakVerticalDetector[:, 1]
                        )
                    )
                    / np.sum(peakVerticalDetector[:, 1]),
                ),
                axis=0,
            )  ##
            uncertaintyFitParamsHD = np.append(
                uncertaintyFitParamsHD,
                np.sqrt(np.diag(covarianceHD))[5 * n : 5 * n + 5],
                axis=0,
            )  ##
            uncertaintyFitParamsVD = np.append(
                uncertaintyFitParamsVD,
                np.sqrt(np.diag(covarianceVD))[5 * n : 5 * n + 5],
                axis=0,
            )  ##
    rawDataLevel1_1.create_dataset(
        "horizontalDetector", dtype="float64", data=patternHorizontalDetector
    )  ## save raw data of the horizontal detector
    rawDataLevel1_1.create_dataset(
        "verticalDetector", dtype="float64", data=patternVerticalDetector
    )  ## save raw data of the vertical detector
    fitLevel1_2["fitParams"].create_dataset(
        "fitParamsHD",
        dtype="float64",
        data=np.reshape(fitParamsHD, (int(np.size(fitParamsHD) / 6), 6)),
    )  ## save parameters of the fit of HD
    fitLevel1_2["fitParams"].create_dataset(
        "fitParamsVD",
        dtype="float64",
        data=np.reshape(fitParamsVD, (int(np.size(fitParamsVD) / 6), 6)),
    )  ## save parameters of the fit of VD
    fitLevel1_2["fitParams"].create_dataset(
        "uncertaintyFitParamsHD",
        dtype="float64",
        data=np.reshape(
            uncertaintyFitParamsHD, (int(np.size(uncertaintyFitParamsHD) / 5), 5)
        ),
    )  ## save uncertainty on the parameters of the fit of HD
    fitLevel1_2["fitParams"].create_dataset(
        "uncertaintyFitParamsVD",
        dtype="float64",
        data=np.reshape(
            uncertaintyFitParamsVD, (int(np.size(uncertaintyFitParamsVD) / 5), 5)
        ),
    )  ## save uncertainty on the parameters of the fit of VD
    calibrantSource = np.loadtxt(
        sourceCalibrantFile
    )  ## open source calibration text file
    curveCalibrationHD[:, 0] = fitLevel1_2["fitParams/fitParamsHD"][:, 1]
    curveCalibrationHD[:, 1] = calibrantSource[: np.sum(nbPeaksInBoxes)]
    ucurveCalibrationHD = fitLevel1_2["fitParams/uncertaintyFitParamsHD"][:, 1]
    fitLevel1_2["curveCalibration"].create_dataset(
        "curveCalibrationHD", dtype="float64", data=curveCalibrationHD
    )  ## curve energy VS channels for horizontal detector
    curveCalibrationVD[:, 0] = fitLevel1_2["fitParams/fitParamsVD"][:, 1]
    curveCalibrationVD[:, 1] = calibrantSource[: np.sum(nbPeaksInBoxes)]
    ucurveCalibrationVD = fitLevel1_2["fitParams/uncertaintyFitParamsVD"][:, 1]
    fitLevel1_2["curveCalibration"].create_dataset(
        "curveCalibrationVD", dtype="float64", data=curveCalibrationVD
    )  ## curve energy VS channels for vertical detector
    calibCoeffsHD, covCalibCoeffsHD = np.polyfit(
        x=curveCalibrationHD[:, 0],
        y=curveCalibrationHD[:, 1],
        deg=2,
        full=False,
        cov="unscaled",
        w=1 / ucurveCalibrationHD,
    )  ## calibration coefficients of the horizontal detector
    calibCoeffsVD, covCalibCoeffsVD = np.polyfit(
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
        data=np.transpose(
            (
                curveCalibrationHD[:, 0],
                np.poly1d(calibCoeffsHD)(curveCalibrationHD[:, 0]),
            )
        ),
    )  ## fitted curve energy VS channels for horizontal detector
    fitLevel1_2["curveCalibration"].create_dataset(
        "fitCurveCalibrationVD",
        dtype="float64",
        data=np.transpose(
            (
                curveCalibrationVD[:, 0],
                np.poly1d(calibCoeffsVD)(curveCalibrationVD[:, 0]),
            )
        ),
    )  ## fitted curve energy VS channels for vertical detector
    fitLevel1_2["curveCalibration"].create_dataset(
        "errorCurveCalibrationHD",
        dtype="float64",
        data=np.transpose(
            (
                curveCalibrationHD[:, 0],
                np.abs(
                    np.poly1d(calibCoeffsHD)(curveCalibrationHD[:, 0])
                    - curveCalibrationHD[:, 1]
                ),
            )
        ),
    )  ## error between fitted and raw curve energy VS channels for horizontal detector
    fitLevel1_2["curveCalibration"].create_dataset(
        "errorCurveCalibrationVD",
        dtype="float64",
        data=np.transpose(
            (
                curveCalibrationVD[:, 0],
                np.abs(
                    np.poly1d(calibCoeffsVD)(curveCalibrationVD[:, 0])
                    - curveCalibrationVD[:, 1]
                ),
            )
        ),
    )  ## error between fitted and raw curve energy VS channels for vertical detector
    # print(f'uncertauntyCalibCoeffsHD = {np.sqrt(np.diag(covCalibCoeffsHD))}')
    # print(f'uncertauntyCalibCoeffsVD = {np.sqrt(np.diag(covCalibCoeffsVD))}')
    fitLevel1_2["calibCoeffs"].create_dataset(
        "calibCoeffsHD", dtype="float64", data=calibCoeffsHD
    )  ## save calibration coefficients of the horizontal detector
    fitLevel1_2["calibCoeffs"].create_dataset(
        "calibCoeffsVD", dtype="float64", data=calibCoeffsVD
    )  ## save calibration coefficients of the Vertical detector
    fitLevel1_2["calibCoeffs"].create_dataset(
        "uncertaintyCalibCoeffsHD",
        dtype="float64",
        data=np.sqrt(np.diag(covCalibCoeffsHD)),
    )  ## save uncertainty on calibration coefficients of the horizontal detector
    fitLevel1_2["calibCoeffs"].create_dataset(
        "uncertaintyCalibCoeffsVD",
        dtype="float64",
        data=np.sqrt(np.diag(covCalibCoeffsVD)),
    )  ## save uncertainty on calibration coefficients of the vertical detector

    h5Save.close()
    return


if __name__ == "__main__":
    run_from_cli(calibEdd)
