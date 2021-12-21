from typing import Sequence
import numpy as np
import h5py
import silx.math.fit
import scipy.optimize
from easistrain.EDD.io import as_nxchar, create_info_group, peak_dataset_data
from easistrain.EDD.utils import (
    process_detector_data,
    run_from_cli,
    splitPseudoVoigt,
    calcBackground,
    guessParameters,
)


def fitEDD(
    fileRead: str,
    fileSave: str,
    sample: str,
    dataset: str,
    scanNumber: int,
    nameHorizontalDetector: str,
    nameVerticalDetector: str,
    positioners: Sequence[str],
    numberOfBoxes: int,
    nbPeaksInBoxes: Sequence[int],
    rangeFitHD: Sequence[int],
    rangeFitVD: Sequence[int],
):
    print(f"Fitting scan n.{scanNumber}")

    with h5py.File(fileRead, "r") as h5Read:  ## Read the h5 file of raw data
        scan_meas = h5Read.get(
            f"{sample}_{dataset}_{scanNumber}.1/measurement",
            default=None,
        )

        if (
            not isinstance(scan_meas, h5py.Group)
            or nameHorizontalDetector not in scan_meas
            or nameVerticalDetector not in scan_meas
        ):
            print("No pattern was saved in this scan")
            return

        h5Save = h5py.File(fileSave, "a")  ## create/append h5 file to save in
        scanGroup = h5Save.create_group(
            f"{sample}_{dataset}_{scanNumber}.1"
        )  ## create the group of the scan wich will contatin all the results of a scan
        positionersGroup = scanGroup.create_group(
            "positioners"
        )  ## positioners subgroup in scan group
        patternHorizontalDetector = h5Read[
            f"{sample}_{dataset}_{scanNumber}.1/measurement/{nameHorizontalDetector}"
        ][
            ()
        ]  ## pattern of horizontal detector
        patternVerticalDetector = h5Read[
            f"{sample}_{dataset}_{scanNumber}.1/measurement/{nameVerticalDetector}"
        ][
            ()
        ]  ## pattern of vertical detector
        twoD_detector_data = (
            np.ndim(patternHorizontalDetector) == 2
            or np.ndim(patternVerticalDetector) == 2
        )

        if twoD_detector_data:
            positionAngles = np.zeros((len(patternHorizontalDetector), 6), "float64")
        else:
            positionAngles = np.zeros((1, 6), "float64")
        for i, positioner in enumerate(positioners):
            pos_data = h5Read[
                f"{sample}_{dataset}_{scanNumber}.1/instrument/positioners/{positioner}"
            ][()]
            positionersGroup.create_dataset(
                positioner,
                dtype="float64",
                data=pos_data,
            )  ## saving all the requested positioners
            if i < 6:
                positionAngles[:, i] = pos_data
            else:
                print("Too many positioners given ! Only 6 are handled for now.")

    rawDataLevel1_1 = scanGroup.create_group(
        "rawData" + "_" + str(dataset) + "_" + str(scanNumber)
    )  ## rawData subgroup in scan group
    fitGroup = scanGroup.create_group("fit")  ## fit subgroup in scan group
    tthPositionsGroup = scanGroup.create_group(
        "tthPositionsGroup"
    )  ## two theta positions subgroup in scan group
    rawDataLevel1_1.create_dataset(
        "horizontalDetector", dtype="float64", data=patternHorizontalDetector
    )  ## save raw data of the horizontal detector
    rawDataLevel1_1.create_dataset(
        "verticalDetector", dtype="float64", data=patternVerticalDetector
    )  ## save raw data of the vertical detector
    if twoD_detector_data:
        for k in range(len(patternHorizontalDetector)):
            fitParamsHD = np.array(())
            fitParamsVD = np.array(())
            uncertaintyFitParamsHD = np.array(())
            uncertaintyFitParamsVD = np.array(())
            pointInScan = fitGroup.create_group(
                f"{str(k).zfill(4)}"
            )  ## create a group of each pattern (point of the scan)
            pointInScan.create_group(
                "fitParams"
            )  ## fit results group for the two detector
            for i in range(numberOfBoxes):
                peakHorizontalDetector = np.transpose(
                    (
                        np.arange(rangeFitHD[2 * i], rangeFitHD[(2 * i) + 1]),
                        patternHorizontalDetector[
                            k, rangeFitHD[2 * i] : rangeFitHD[(2 * i) + 1]
                        ],
                    )
                )  ## peak of the horizontal detector
                peakVerticalDetector = np.transpose(
                    (
                        np.arange(rangeFitVD[2 * i], rangeFitVD[(2 * i) + 1]),
                        patternVerticalDetector[
                            k, rangeFitVD[2 * i] : rangeFitVD[(2 * i) + 1]
                        ],
                    )
                )  ## peak of the vertical detector
                backgroundHorizontalDetector = silx.math.fit.strip(
                    data=peakHorizontalDetector[:, 1],
                    w=5,
                    niterations=4000,
                    factor=1,
                    anchors=None,
                )  ## stripped background of the horizontal detector (obtained by stripping the yData)
                backgroundVerticalDetector = silx.math.fit.strip(
                    data=peakVerticalDetector[:, 1],
                    w=5,
                    niterations=4000,
                    factor=1,
                    anchors=None,
                )  ## stripped background of the vertical detector (obtained by stripping the yData)
                fitLine = pointInScan.create_group(
                    f"fitLine_{str(i).zfill(4)}"
                )  ## create group for each range of peak(s)
                fitLine.create_dataset(
                    "rawHorizontalDetector",
                    dtype="float64",
                    data=peakHorizontalDetector,
                )  ## create dataset for raw data of each calibration peak
                fitLine.create_dataset(
                    "rawVerticalDetector", dtype="f", data=peakVerticalDetector
                )  ## create dataset for raw data of each calibration peak
                peaksGuessHD, peaksIndexHD = guessParameters(
                    peakHorizontalDetector[:, 0],
                    peakHorizontalDetector[:, 1] - backgroundHorizontalDetector,
                    nbPeaksInBoxes[i],
                    withBounds=True,
                )  ## guess fit parameters for HD
                peaksGuessVD, peaksIndexVD = guessParameters(
                    peakVerticalDetector[:, 0],
                    peakVerticalDetector[:, 1] - backgroundVerticalDetector,
                    nbPeaksInBoxes[i],
                    withBounds=True,
                )  ## guess fit parameters for VD
                yCalculatedBackgroundHD = calcBackground(
                    peakHorizontalDetector[:, 0],
                    peakHorizontalDetector[:, 1],
                    peaksGuessHD[-1],
                    peaksGuessHD[2],
                    peaksIndexHD,
                )  ## calculated ybackground of the horizontal detector
                fitLine.create_dataset(
                    "backgroundHorizontalDetector",
                    dtype="float64",
                    data=np.transpose(
                        (peakHorizontalDetector[:, 0], yCalculatedBackgroundHD)
                    ),
                )  ## create dataset for background of each calibration peak for HD
                fitLine.create_dataset(
                    "bgdSubsDataHorizontalDetector",
                    dtype="float64",
                    data=np.transpose(
                        (
                            peakHorizontalDetector[:, 0],
                            peakHorizontalDetector[:, 1] - yCalculatedBackgroundHD,
                        )
                    ),
                )  ## create dataset for HD raw data after subst of background
                initialGuessHD = np.zeros(5 * nbPeaksInBoxes[i])
                minBoundsHD = np.array(())
                maxBoundsHD = np.array(())
                for n in range(nbPeaksInBoxes[i]):
                    initialGuessHD[5 * n] = peaksGuessHD[3 * n]
                    initialGuessHD[5 * n + 1] = peaksGuessHD[3 * n + 1]
                    initialGuessHD[5 * n + 2] = peaksGuessHD[3 * n + 2]
                    initialGuessHD[5 * n + 3] = peaksGuessHD[3 * n + 2]
                    initialGuessHD[5 * n + 4] = 0.5
                    appendMinBoundsHD = np.array(
                        (
                            [
                                np.amin(peakHorizontalDetector[:, 1]),
                                initialGuessHD[5 * n + 1]
                                - 3 * initialGuessHD[5 * n + 2],
                                0,
                                0,
                                0,
                            ]
                        )
                    )  # minimum bounds of the parametrs solution (H, C, FWHM1, FWHM2, eta) to append for the horizontal detector
                    appendMaxBoundsHD = np.array(
                        (
                            [
                                np.amax(peakHorizontalDetector[:, 1]),
                                initialGuessHD[5 * n + 1]
                                + 3 * initialGuessHD[5 * n + 2],
                                len(peakHorizontalDetector[:, 0]) / 2,
                                len(peakHorizontalDetector[:, 0]) / 2,
                                1,
                            ]
                        )
                    )  # maximum bounds of the parametrs solution (H, C, FWHM1, FWHM2, eta)to append for the horizontal detector
                    minBoundsHD = np.append(
                        minBoundsHD, appendMinBoundsHD
                    )  # minimum bounds of the parametrs solution (H, C, FWHM1, FWHM2, eta) for the horizontal detector
                    maxBoundsHD = np.append(
                        maxBoundsHD, appendMaxBoundsHD
                    )  # maximum bounds of the parametrs solution (H, C, FWHM1, FWHM2, eta) for the horizontal detector
                optimal_parametersHD, covarianceHD = scipy.optimize.curve_fit(
                    f=splitPseudoVoigt,
                    xdata=peakHorizontalDetector[:, 0],
                    ydata=peakHorizontalDetector[:, 1] - yCalculatedBackgroundHD,
                    p0=initialGuessHD,
                    sigma=None,
                    bounds=(minBoundsHD, maxBoundsHD),
                    maxfev=10000,
                )  ## fit of the peak of the Horizontal detector
                fitLine.create_dataset(
                    "fitHorizontalDetector",
                    dtype="float64",
                    data=np.transpose(
                        (
                            peakHorizontalDetector[:, 0],
                            splitPseudoVoigt(
                                peakHorizontalDetector[:, 0], optimal_parametersHD
                            )
                            + yCalculatedBackgroundHD,
                        )
                    ),
                )  ## fitted data of the horizontal detector
                fitLine.create_dataset(
                    "errorHorizontalDetector",
                    dtype="float64",
                    data=np.transpose(
                        (
                            peakHorizontalDetector[:, 0],
                            np.absolute(
                                splitPseudoVoigt(
                                    peakHorizontalDetector[:, 0],
                                    optimal_parametersHD,
                                )
                                + yCalculatedBackgroundHD
                                - peakHorizontalDetector[:, 1]
                            ),
                        )
                    ),
                )  ## error of the horizontal detector
                for n in range(nbPeaksInBoxes[i]):
                    fitParamsHD = np.append(
                        fitParamsHD,
                        np.append(
                            optimal_parametersHD[5 * n : 5 * n + 5],
                            100
                            * np.sum(
                                np.absolute(
                                    splitPseudoVoigt(
                                        peakHorizontalDetector[:, 0],
                                        optimal_parametersHD,
                                    )
                                    + backgroundHorizontalDetector
                                    - peakHorizontalDetector[:, 1]
                                )
                            )
                            / np.sum(peakHorizontalDetector[:, 1]),
                        ),
                        axis=0,
                    )  ##
                    uncertaintyFitParamsHD = np.append(
                        uncertaintyFitParamsHD,
                        np.sqrt(np.diag(covarianceHD))[5 * n : 5 * n + 5],
                        axis=0,
                    )  ##
                yCalculatedBackgroundVD = calcBackground(
                    peakVerticalDetector[:, 0],
                    peakVerticalDetector[:, 1],
                    peaksGuessVD[-1],
                    peaksGuessVD[2],
                    peaksIndexVD,
                )  ## calculated ybackground of the vertical detector
                fitLine.create_dataset(
                    "backgroundVerticalDetector",
                    dtype="float64",
                    data=np.transpose(
                        (peakVerticalDetector[:, 0], yCalculatedBackgroundVD)
                    ),
                )  ## create dataset for background of each calibration peak for VD
                fitLine.create_dataset(
                    "bgdSubsDataVerticalDetector",
                    dtype="float64",
                    data=np.transpose(
                        (
                            peakVerticalDetector[:, 0],
                            peakVerticalDetector[:, 1] - yCalculatedBackgroundVD,
                        )
                    ),
                )  ## create dataset for VD raw data after subst of background
                initialGuessVD = np.zeros(5 * nbPeaksInBoxes[i])
                minBoundsVD = np.array(())
                maxBoundsVD = np.array(())
                for n in range(nbPeaksInBoxes[i]):
                    initialGuessVD[5 * n] = peaksGuessVD[3 * n]
                    initialGuessVD[5 * n + 1] = peaksGuessVD[3 * n + 1]
                    initialGuessVD[5 * n + 2] = peaksGuessVD[3 * n + 2]
                    initialGuessVD[5 * n + 3] = peaksGuessVD[3 * n + 2]
                    initialGuessVD[5 * n + 4] = 0.5
                    appendMinBoundsVD = np.array(
                        (
                            [
                                np.amin(peakVerticalDetector[:, 1]),
                                initialGuessVD[5 * n + 1]
                                - 3 * initialGuessVD[5 * n + 2],
                                0,
                                0,
                                0,
                            ]
                        )
                    )  # minimum bounds of the parametrs solution (H, C, FWHM1, FWHM2, eta) to append for the vertical detector
                    appendMaxBoundsVD = np.array(
                        (
                            [
                                np.amax(peakVerticalDetector[:, 1]),
                                initialGuessVD[5 * n + 1]
                                + 3 * initialGuessVD[5 * n + 2],
                                len(peakVerticalDetector[:, 0]) / 2,
                                len(peakVerticalDetector[:, 0]) / 2,
                                1,
                            ]
                        )
                    )  # maximum bounds of the parametrs solution (H, C, FWHM1, FWHM2, eta)to append for the vertical detector
                    minBoundsVD = np.append(
                        minBoundsVD, appendMinBoundsVD
                    )  # minimum bounds of the parametrs solution (H, C, FWHM1, FWHM2, eta) for the vertical detector
                    maxBoundsVD = np.append(
                        maxBoundsVD, appendMaxBoundsVD
                    )  # maximum bounds of the parametrs solution (H, C, FWHM1, FWHM2, eta) for the vertical detector
                optimal_parametersVD, covarianceVD = scipy.optimize.curve_fit(
                    f=splitPseudoVoigt,
                    xdata=peakVerticalDetector[:, 0],
                    ydata=peakVerticalDetector[:, 1] - yCalculatedBackgroundVD,
                    p0=initialGuessVD,
                    sigma=None,
                    bounds=(minBoundsVD, maxBoundsVD),
                    maxfev=10000,
                )  ## fit of the peak of the Vertical detector
                fitLine.create_dataset(
                    "fitVerticalDetector",
                    dtype="float64",
                    data=np.transpose(
                        (
                            peakVerticalDetector[:, 0],
                            splitPseudoVoigt(
                                peakVerticalDetector[:, 0], optimal_parametersVD
                            )
                            + yCalculatedBackgroundVD,
                        )
                    ),
                )  ## fitted data of the vertical detector
                fitLine.create_dataset(
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
                for n in range(nbPeaksInBoxes[i]):
                    fitParamsVD = np.append(
                        fitParamsVD,
                        np.append(
                            optimal_parametersVD[5 * n : 5 * n + 5],
                            100
                            * np.sum(
                                np.absolute(
                                    splitPseudoVoigt(
                                        peakVerticalDetector[:, 0],
                                        optimal_parametersVD,
                                    )
                                    + backgroundVerticalDetector
                                    - peakVerticalDetector[:, 1]
                                )
                            )
                            / np.sum(peakVerticalDetector[:, 1]),
                        ),
                        axis=0,
                    )  ##
                    uncertaintyFitParamsVD = np.append(
                        uncertaintyFitParamsVD,
                        np.sqrt(np.diag(covarianceVD))[5 * n : 5 * n + 5],
                        axis=0,
                    )  ##
            pointInScan["fitParams"].create_dataset(
                "fitParamsHD",
                dtype="float64",
                data=np.reshape(fitParamsHD, (int(np.size(fitParamsHD) / 6), 6)),
            )  ## save parameters of the fit of HD
            pointInScan["fitParams"].create_dataset(
                "uncertaintyFitParamsHD",
                dtype="float64",
                data=np.reshape(
                    uncertaintyFitParamsHD,
                    (int(np.size(uncertaintyFitParamsHD) / 5), 5),
                ),
            )  ## save uncertainty on the parameters of the fit of HD

            pointInScan["fitParams"].create_dataset(
                "fitParamsVD",
                dtype="float64",
                data=np.reshape(fitParamsVD, (int(np.size(fitParamsVD) / 6), 6)),
            )  ## save parameters of the fit of VD
            pointInScan["fitParams"].create_dataset(
                "uncertaintyFitParamsVD",
                dtype="float64",
                data=np.reshape(
                    uncertaintyFitParamsVD,
                    (int(np.size(uncertaintyFitParamsVD) / 5), 5),
                ),
            )  ## save uncertainty on the parameters of the fit of VD
            for peakNumber in range(np.sum(nbPeaksInBoxes)):
                if f"peak_{str(peakNumber).zfill(4)}" not in tthPositionsGroup.keys():
                    peakDataset = tthPositionsGroup.create_dataset(
                        f"peak_{str(peakNumber).zfill(4)}",
                        dtype="float64",
                        data=np.zeros(
                            (2 * len(patternHorizontalDetector), 13), "float64"
                        ),
                    )  ## create a group for each peak in tthPositionGroup
                    uncertaintyPeakDataset = tthPositionsGroup.create_dataset(
                        f"uncertaintyPeak_{str(peakNumber).zfill(4)}",
                        dtype="float64",
                        data=np.zeros(
                            (2 * len(patternHorizontalDetector), 13), "float64"
                        ),
                    )  ## create a group for each peak in tthPositionGroup
                else:
                    peakDataset = tthPositionsGroup[f"peak_{str(peakNumber).zfill(4)}"]
                    uncertaintyPeakDataset = tthPositionsGroup[
                        f"uncertaintyPeak_{str(peakNumber).zfill(4)}"
                    ]
                peakDataset[2 * k, 0:6] = positionAngles[
                    k, 0:6
                ]  ## coordinates of the point in the insrument reference and phi, chi and omega angles of the instrument
                uncertaintyPeakDataset[2 * k, 0:6] = positionAngles[
                    k, 0:6
                ]  ## coordinates of the point in the insrument reference and phi, chi and omega angles of the instrument (for the moment I set it to zero) in the uncertainty dataset
                peakDataset[
                    2 * k, 6
                ] = (
                    -90
                )  ## delta angle of the horizontal detector (debye scherer ring angle)
                uncertaintyPeakDataset[
                    2 * k, 6
                ] = (
                    -90
                )  ## delta angle of the horizontal detector (debye scherer ring angle) (I set the uncertainty to zero for the moment)
                peakDataset[
                    2 * k, 7
                ] = 0  ## theta angle (diffraction fixed angle) of the horizontal detector (I suppose that it is zero as we work at high energy and the angle is fixed to 2.5 deg)
                uncertaintyPeakDataset[
                    2 * k, 7
                ] = 0  ## uncertainty of the theta angle (diffraction fixed angle) of the horizontal detector (I suppose that it is zero as we work at high energy and the angle is fixed to 2.5 deg)
                peakDataset[2 * k, 8] = pointInScan["fitParams/fitParamsHD"][
                    (peakNumber, 1)
                ]  ## peak position of HD
                uncertaintyPeakDataset[2 * k, 8] = pointInScan[
                    "fitParams/uncertaintyFitParamsHD"
                ][
                    (peakNumber, 1)
                ]  ## uncertainty of the peak position of HD
                peakDataset[2 * k, 9] = pointInScan["fitParams/fitParamsHD"][
                    (peakNumber, 0)
                ]  ## peak intensity of HD (maximum intensity)
                uncertaintyPeakDataset[2 * k, 9] = pointInScan[
                    "fitParams/uncertaintyFitParamsHD"
                ][
                    (peakNumber, 0)
                ]  ## uncertainty of the peak intensity of HD (maximum intensity)
                peakDataset[2 * k, 10] = 0.5 * (
                    pointInScan["fitParams/fitParamsHD"][(peakNumber, 2)]
                    + pointInScan["fitParams/fitParamsHD"][(peakNumber, 3)]
                )  ## FWHM of HD (mean of the FWHM at left and right as I am using an assymetric function)
                uncertaintyPeakDataset[2 * k, 10] = 0.5 * (
                    pointInScan["fitParams/uncertaintyFitParamsHD"][(peakNumber, 2)]
                    + pointInScan["fitParams/uncertaintyFitParamsHD"][(peakNumber, 3)]
                )  ## uncertainty if the FWHM of HD (mean of the FWHM at left and right as I am using an assymetric function)
                peakDataset[2 * k, 11] = pointInScan["fitParams/fitParamsHD"][
                    (peakNumber, 4)
                ]  ## shape factor of HD (contribution of the lorentzian function)
                uncertaintyPeakDataset[2 * k, 11] = pointInScan[
                    "fitParams/uncertaintyFitParamsHD"
                ][
                    (peakNumber, 4)
                ]  ## uncertainty of the shape factor of HD (contribution of the lorentzian function)
                peakDataset[2 * k, 12] = pointInScan["fitParams/fitParamsHD"][
                    (peakNumber, 5)
                ]  ## Rw factor of HD (goodness of the fit)
                uncertaintyPeakDataset[
                    2 * k, 12
                ] = 0  ## set it to zero no utility of this thing /uncertainty of Rw factor of HD (goodness of the fit)
                peakDataset[2 * k + 1, 0:6] = positionAngles[
                    k, 0:6
                ]  ## coordinates of the point in the insrument reference and phi, chi and omega angles of the instrument
                uncertaintyPeakDataset[2 * k + 1, 0:6] = positionAngles[
                    k, 0:6
                ]  ## uncertaiinty on the coordinates of the point in the insrument reference and phi, chi and omega angles of the instrument (I set it to zero for the moment)
                peakDataset[
                    2 * k + 1, 6
                ] = 0  ## delta angle of the vertical detector (debye scherer ring angle)
                uncertaintyPeakDataset[
                    2 * k + 1, 6
                ] = 0  ## uncertainty of the delta angle of the vertical detector (debye scherer ring angle) (I set the uncertainty to zero for the moment)
                peakDataset[
                    2 * k + 1, 7
                ] = 0  ## theta angle (diffraction fixed angle) of the vertical detector (I suppose that it is zero as we work a small angle fixed to 2.5 deg)
                uncertaintyPeakDataset[
                    2 * k + 1, 7
                ] = 0  ## uncertainty of the theta angle (diffraction fixed angle) of the vertical detector (I suppose that it is zero as we work with a small angle fixed to 2.5 deg)
                peakDataset[2 * k + 1, 8] = pointInScan["fitParams/fitParamsVD"][
                    (peakNumber, 1)
                ]  ## peak position of VD
                uncertaintyPeakDataset[2 * k + 1, 8] = pointInScan[
                    "fitParams/uncertaintyFitParamsVD"
                ][
                    (peakNumber, 1)
                ]  ## uncertainty of the peak position of VD
                peakDataset[2 * k + 1, 9] = pointInScan["fitParams/fitParamsVD"][
                    (peakNumber, 0)
                ]  ## peak intensity of VD (maximum intensity)
                uncertaintyPeakDataset[2 * k + 1, 9] = pointInScan[
                    "fitParams/uncertaintyFitParamsVD"
                ][
                    (peakNumber, 0)
                ]  ## uncertainty of the peak intensity of VD (maximum intensity)
                peakDataset[2 * k + 1, 10] = 0.5 * (
                    pointInScan["fitParams/fitParamsVD"][(peakNumber, 2)]
                    + pointInScan["fitParams/fitParamsVD"][(peakNumber, 3)]
                )  ## FWHM of VD (mean of the FWHM at left and right as I am using an assymetric function)
                uncertaintyPeakDataset[2 * k + 1, 10] = 0.5 * (
                    pointInScan["fitParams/uncertaintyFitParamsVD"][(peakNumber, 2)]
                    + pointInScan["fitParams/uncertaintyFitParamsVD"][(peakNumber, 3)]
                )  ## uncertainty of the FWHM of VD (mean of the FWHM at left and right as I am using an assymetric function)
                peakDataset[2 * k + 1, 11] = pointInScan["fitParams/fitParamsVD"][
                    (peakNumber, 4)
                ]  ## shape factor of VD (contribution of the lorentzian function)
                uncertaintyPeakDataset[2 * k + 1, 11] = pointInScan[
                    "fitParams/uncertaintyFitParamsVD"
                ][
                    (peakNumber, 4)
                ]  ## uncertainty of the shape factor of VD (contribution of the lorentzian function)
                peakDataset[2 * k + 1, 12] = pointInScan["fitParams/fitParamsVD"][
                    (peakNumber, 5)
                ]  ## Rw factor of VD (goodness of the fit)
                uncertaintyPeakDataset[
                    2 * k + 1, 12
                ] = 0  ## No utility of this thing, I set it to zero/ uncertainty of the Rw factor of VD (goodness of the fit)
            if "infoPeak" not in tthPositionsGroup.keys():
                tthPositionsGroup.create_dataset(
                    "infoPeak",
                    dtype=h5py.string_dtype(encoding="utf-8"),
                    data=f"{positioners}, delta, thetha, position in channel, Intenstity, FWHM, shape factor, goodness factor",
                )  ## create info about dataset saved for each peak in tthPositionGroup
        # print(positionAngles)
    else:
        fitParams = {"horizontal": np.array(()), "vertical": np.array(())}
        uncertaintyFitParams = {"horizontal": np.array(()), "vertical": np.array(())}
        pointInScan = fitGroup.create_group(
            f"{str(0).zfill(4)}"
        )  ## create a group of each pattern (point of the scan)
        fitParamsGroup = pointInScan.create_group(
            "fitParams"
        )  ## fit results group for the two detector
        for i, nb_peaks in enumerate(nbPeaksInBoxes):
            fitLine = pointInScan.create_group(
                f"fitLine_{str(i).zfill(4)}"
            )  ## create group for each range of peak(s)

            for detector in ["horizontal", "vertical"]:
                range_fit = (
                    rangeFitHD if detector == "horizontal" else rangeFitVD
                )  # To be improved
                pattern = (
                    patternHorizontalDetector
                    if detector == "horizontal"
                    else patternVerticalDetector
                )  # To be improved
                assert isinstance(pattern, np.ndarray)
                (
                    peakDetector,
                    yCalculatedBackground,
                    optimal_parameters,
                    boxFitParams,
                    uncertaintyBoxFitParams,
                ) = process_detector_data(
                    fit_min=range_fit[2 * i],
                    fit_max=range_fit[(2 * i) + 1],
                    input_data=pattern,
                    nb_peaks=nb_peaks,
                )
                horizontal_channels = peakDetector[:, 0]
                raw_horizontal_data = peakDetector[:, 1]
                horizontal_fit = (
                    splitPseudoVoigt(horizontal_channels, optimal_parameters)
                    + yCalculatedBackground
                )

                detectorGroup = fitLine.create_group(detector)
                detectorGroup.create_dataset(
                    "channels",
                    dtype="float64",
                    data=horizontal_channels,
                )
                detectorGroup.create_dataset(
                    "raw_data",
                    dtype="float64",
                    data=raw_horizontal_data,
                )
                detectorGroup.create_dataset(
                    "background",
                    dtype="float64",
                    data=yCalculatedBackground,
                )
                detectorGroup.create_dataset(
                    "data - background",
                    dtype="float64",
                    data=raw_horizontal_data - yCalculatedBackground,
                )
                detectorGroup.create_dataset(
                    "fitted_data",
                    dtype="float64",
                    data=horizontal_fit,
                )
                detectorGroup.create_dataset(
                    "error",
                    dtype="float64",
                    data=(horizontal_fit - raw_horizontal_data),
                )  ## error of the horizontal detector

                # NeXus
                detectorGroup.attrs["NX_class"] = as_nxchar("NXdata")
                detectorGroup.attrs["auxiliary_signals"] = as_nxchar(["raw_data"])
                detectorGroup.attrs["errors"] = as_nxchar("error")
                detectorGroup.attrs["signal"] = as_nxchar("fitted_data")
                detectorGroup.attrs["axes"] = as_nxchar("channels")

                # Accumulate fit parameters of this box
                fitParams[detector] = np.append(fitParams[detector], boxFitParams)
                uncertaintyFitParams[detector] = np.append(
                    uncertaintyFitParams[detector], uncertaintyBoxFitParams
                )

        # End of detector processing
        savedFitParamsHD = np.reshape(
            fitParams["horizontal"], (int(np.size(fitParams["horizontal"]) / 6), 6)
        )
        fitParamsGroup.create_dataset(
            "fitParamsHD",
            dtype="float64",
            data=savedFitParamsHD,
        )  ## save parameters of the fit of HD
        savedUncertaintyFitParamsHD = np.reshape(
            uncertaintyFitParams["horizontal"],
            (int(np.size(uncertaintyFitParams["horizontal"]) / 5), 5),
        )
        fitParamsGroup.create_dataset(
            "uncertaintyFitParamsHD",
            dtype="float64",
            data=savedUncertaintyFitParamsHD,
        )  ## save uncertainty on the parameters of the fit of HD

        savedFitParamsVD = np.reshape(
            fitParams["vertical"], (int(np.size(fitParams["vertical"]) / 6), 6)
        )
        fitParamsGroup.create_dataset(
            "fitParamsVD",
            dtype="float64",
            data=savedFitParamsVD,
        )  ## save parameters of the fit of VD
        savedUncertaintyFitParamsVD = np.reshape(
            uncertaintyFitParams["vertical"],
            (int(np.size(uncertaintyFitParams["vertical"]) / 5), 5),
        )
        fitParamsGroup.create_dataset(
            "uncertaintyFitParamsVD",
            dtype="float64",
            data=savedUncertaintyFitParamsVD,
        )  ## save uncertainty on the parameters of the fit of VD
        for peakNumber in range(np.sum(nbPeaksInBoxes)):
            if f"peak_{str(peakNumber).zfill(4)}" not in tthPositionsGroup.keys():
                peakDataset = tthPositionsGroup.create_dataset(
                    f"peak_{str(peakNumber).zfill(4)}",
                    dtype="float64",
                    data=np.zeros((2, 13), "float64"),
                )  ## create a dataset for each peak in tthPositionGroup
                uncertaintyPeakDataset = tthPositionsGroup.create_dataset(
                    f"uncertaintyPeak_{str(peakNumber).zfill(4)}",
                    dtype="float64",
                    data=np.zeros((2, 13), "float64"),
                )  ## create a dataset for uncertainty for each peak in tthPositionGroup
            else:
                peakDataset = tthPositionsGroup[f"peak_{str(peakNumber).zfill(4)}"]
                assert isinstance(peakDataset, h5py.Dataset)
                uncertaintyPeakDataset = tthPositionsGroup[
                    f"uncertaintyPeak_{str(peakNumber).zfill(4)}"
                ]
                assert isinstance(uncertaintyPeakDataset, h5py.Dataset)
            peakDataset[0] = peak_dataset_data(
                positionAngles, savedFitParamsHD[peakNumber]
            )
            peakDataset[1] = peak_dataset_data(
                positionAngles, savedFitParamsVD[peakNumber]
            )
            uncertaintyPeakDataset[0] = peak_dataset_data(
                positionAngles, savedUncertaintyFitParamsHD[peakNumber]
            )
            uncertaintyPeakDataset[1] = peak_dataset_data(
                positionAngles, savedUncertaintyFitParamsVD[peakNumber]
            )
        if "infoPeak" not in tthPositionsGroup.keys():
            tthPositionsGroup.create_dataset(
                "infoPeak",
                dtype=h5py.string_dtype(encoding="utf-8"),
                data=f"{positioners}, delta, theta, position in channel, Intenstity, FWHM, shape factor, goodness factor",
            )  ## create info about dataset saved for each peak in tthPositionGroup

    create_info_group(
        scanGroup,
        fileRead,
        fileSave,
        sample,
        dataset,
        scanNumber,
        nameHorizontalDetector,
        nameVerticalDetector,
        numberOfBoxes,
        nbPeaksInBoxes,
        rangeFitHD,
        rangeFitVD,
        positioners,
    )

    h5Save.close()

    return


def fitEDD_with_scan_number_parse(**config):
    """Wrapper function to allow scanNumber to be a list or a slice."""
    n_scan_arg = config.pop("scanNumber")
    if isinstance(n_scan_arg, int):
        fitEDD(**config, scanNumber=n_scan_arg)
    elif isinstance(n_scan_arg, list):
        for i in n_scan_arg:
            fitEDD_with_scan_number_parse(**config, scanNumber=i)
    elif isinstance(n_scan_arg, str):
        if ":" in n_scan_arg:
            min_scan, max_scan = n_scan_arg.split(":")
            for i in range(int(min_scan), int(max_scan)):
                fitEDD(**config, scanNumber=i)
        else:
            fitEDD(**config, scanNumber=int(n_scan_arg))
    else:
        raise ValueError(f"Unrecognized value for scanNumber: {n_scan_arg}")


if __name__ == "__main__":
    run_from_cli(fitEDD_with_scan_number_parse)
