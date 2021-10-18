import h5py
import numpy as np
import silx.math.fit
import silx.math.fit.peaks


# fileRead = '/home/esrf/slim/data/ihme10/id15/TiC_Calib/ihme10_TiC_calib.h5'
# filesave = '/home/esrf/slim/easistrain/easistrain/EDD/Results_ihme10_TiC_calib.h5'
# sample = 'TiC_calib'
# dataset = '0001'
# scanNumber = '4'
# horizontalDetector = 'mca2_det0'
# verticalDetector = 'mca2_det1'
# numberOfPeaks = 1
# rangeFit = [680,820]
# doublet = [1]


def splitPseudoVoigt(xData, *params):
    return silx.math.fit.sum_splitpvoigt(xData, *params)


def guessParameters(yData, counterOfPeak, doublet):
    fwhmGuess = silx.math.fit.peaks.guess_fwhm(yData)
    peaksGuess = silx.math.fit.peaks.peak_search(
        yData,
        fwhmGuess,
        sensitivity=1,
        begin_index=None,
        end_index=None,
        debug=False,
        relevance_info=False,
    )  ## index of the peak
    if np.size(peaksGuess) > doublet[counterOfPeak]:
        # print(peaksGuess[np.argsort(yData[peaksGuess[:].astype(int)])])
        peaksGuess = peaksGuess[np.argsort(yData[peaksGuess[:].astype(int)])][
            -doublet[counterOfPeak] :
        ]
    print(peaksGuess)
    return fwhmGuess, peaksGuess


def angleCalibrationEDD(
    fileRead,
    fileSave,
    sample,
    dataset,
    scanNumber,
    horizontalDetector,
    verticalDetector,
    numberOfPeaks,
    doublet,
    rangeFit,
):
    with h5py.File(fileRead, "r") as h5Read:  ## Read the h5 file of raw data
        patternHorizontalDetector = h5Read[
            sample
            + "_"
            + str(dataset)
            + "_"
            + str(scanNumber)
            + ".1/measurement/"
            + horizontalDetector
        ][
            ()
        ]  ## pattern of horizontal detector
        patternVerticalDetector = h5Read[
            sample
            + "_"
            + str(dataset)
            + "_"
            + str(scanNumber)
            + ".1/measurement/"
            + verticalDetector
        ][
            ()
        ]  ## pattern of vertical detector

    h5Save = h5py.File(fileSave, "a")  ## create/append h5 file to save in
    if not "angleCalibration" in h5Save.keys():
        angleCalibrationLevel1 = h5Save.create_group(
            "angleCalibration"
        )  ## angleCalibration group
    else:
        angleCalibrationLevel1 = h5Save["angleCalibration"]
    rawDataLevel1_1 = angleCalibrationLevel1.create_group(
        "rawData" + "_" + str(dataset) + "_" + str(scanNumber)
    )  ## rawData subgroup in calibration group
    fitLevel1_2 = angleCalibrationLevel1.create_group(
        "fit" + "_" + str(dataset) + "_" + str(scanNumber)
    )  ## fit subgroup in calibration group
    fitLevel1_2.create_group("fitParams")  ## fit results group for the two detector
    fitParamsHD = np.array(())
    fitParamsVD = np.array(())
    uncertaintyFitParamsHD = np.array(())
    uncertaintyFitParamsVD = np.array(())
    for i in range(numberOfPeaks):
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
            niterations=4000,
            factor=1,
            anchors=None,
        )  ## background of the horizontal detector
        backgroundVerticalDetector = silx.math.fit.strip(
            data=peakVerticalDetector[:, 1],
            w=5,
            niterations=4000,
            factor=1,
            anchors=None,
        )  ## background of the vertical detector
        fitLevel1_2.create_group(
            f"fitLine_{str(i)}"
        )  ## create group for each calibration peak
        fitLevel1_2[f"fitLine_{str(i)}"].create_dataset(
            "rawHorizontalDetector", dtype="f", data=peakHorizontalDetector
        )  ## create dataset for raw data of each calibration peak
        fitLevel1_2[f"fitLine_{str(i)}"].create_dataset(
            "rawVerticalDetector", dtype="f", data=peakVerticalDetector
        )  ## create dataset for raw data of each calibration peak
        fitLevel1_2[f"fitLine_{str(i)}"].create_dataset(
            "backgroundHorizontalDetector",
            dtype="f",
            data=np.transpose(
                (peakHorizontalDetector[:, 0], backgroundHorizontalDetector)
            ),
        )  ## create dataset for background of each calibration peak
        fitLevel1_2[f"fitLine_{str(i)}"].create_dataset(
            "backgroundVerticalDetector",
            dtype="f",
            data=np.transpose((peakVerticalDetector[:, 0], backgroundVerticalDetector)),
        )  ## create dataset for background of each calibration peak
        fitLevel1_2[f"fitLine_{str(i)}"].create_dataset(
            "bgdSubsDataHorizontalDetector",
            dtype="f",
            data=np.transpose(
                (
                    peakHorizontalDetector[:, 0],
                    peakHorizontalDetector[:, 1] - backgroundHorizontalDetector,
                )
            ),
        )  ## create dataset for HD raw data after subst of background
        fitLevel1_2[f"fitLine_{str(i)}"].create_dataset(
            "bgdSubsDataVerticalDetector",
            dtype="f",
            data=np.transpose(
                (
                    peakVerticalDetector[:, 0],
                    peakVerticalDetector[:, 1] - backgroundVerticalDetector,
                )
            ),
        )  ## create dataset for VD raw data after subst of background
        fwhmGuessHD, peaksGuessHD = guessParameters(
            peakHorizontalDetector[:, 1] - backgroundHorizontalDetector, i, doublet
        )  ## guess fit parameters for HD
        fwhmGuessVD, peaksGuessVD = guessParameters(
            peakVerticalDetector[:, 1] - backgroundVerticalDetector, i, doublet
        )  ## guess fit parameters for VD
        initialGuessHD = np.zeros(5 * doublet[i])
        initialGuessVD = np.zeros(5 * doublet[i])
        for n in range(doublet[i]):
            initialGuessHD[5 * n] = peakHorizontalDetector[:, 1][int(peaksGuessHD[n])]
            -backgroundHorizontalDetector[int(peaksGuessHD[n])]
            initialGuessHD[5 * n + 1] = peakHorizontalDetector[:, 0][
                int(peaksGuessHD[n])
            ]
            initialGuessHD[5 * n + 2] = fwhmGuessHD
            initialGuessHD[5 * n + 3] = fwhmGuessHD
            initialGuessHD[5 * n + 4] = 0.5
            initialGuessVD[5 * n] = peakVerticalDetector[:, 1][int(peaksGuessVD[n])]
            -backgroundVerticalDetector[int(peaksGuessVD[n])]
            initialGuessVD[5 * n + 1] = peakVerticalDetector[:, 0][int(peaksGuessVD[n])]
            initialGuessVD[5 * n + 2] = fwhmGuessVD
            initialGuessVD[5 * n + 3] = fwhmGuessVD
            initialGuessVD[5 * n + 4] = 0.5
        optimal_parametersHD, covarianceHD, infodictHD = silx.math.fit.leastsq(
            model=splitPseudoVoigt,
            xdata=peakHorizontalDetector[:, 0],
            ydata=peakHorizontalDetector[:, 1] - backgroundHorizontalDetector,
            p0=initialGuessHD,
            sigma=np.sqrt(
                np.abs(peakHorizontalDetector[:, 1] - backgroundHorizontalDetector) + 1
            ),
            full_output=True,
            max_iter=1000,
        )  ## fit of the peak of the Horizontal detector
        optimal_parametersVD, covarianceVD, infodictVD = silx.math.fit.leastsq(
            model=splitPseudoVoigt,
            xdata=peakVerticalDetector[:, 0],
            ydata=peakVerticalDetector[:, 1] - backgroundVerticalDetector,
            p0=initialGuessVD,
            sigma=np.sqrt(
                np.abs(peakVerticalDetector[:, 1] - backgroundVerticalDetector) + 1
            ),
            full_output=True,
            max_iter=1000,
        )  ## fit of the peak of the Vertical detector
        fitLevel1_2[f"fitLine_{str(i)}"].create_dataset(
            "fitHorizontalDetector",
            dtype="f",
            data=np.transpose(
                (
                    peakHorizontalDetector[:, 0],
                    splitPseudoVoigt(peakHorizontalDetector[:, 0], optimal_parametersHD)
                    + backgroundHorizontalDetector,
                )
            ),
        )  ## fitted data of the horizontal detector
        fitLevel1_2[f"fitLine_{str(i)}"].create_dataset(
            "fitVerticalDetector",
            dtype="f",
            data=np.transpose(
                (
                    peakVerticalDetector[:, 0],
                    splitPseudoVoigt(peakVerticalDetector[:, 0], optimal_parametersVD)
                    + backgroundVerticalDetector,
                )
            ),
        )  ## fitted data of the vertical detector
        fitLevel1_2[f"fitLine_{str(i)}"].create_dataset(
            "errorHorizontalDetector",
            dtype="f",
            data=np.transpose(
                (
                    peakHorizontalDetector[:, 0],
                    np.absolute(
                        splitPseudoVoigt(
                            peakHorizontalDetector[:, 0], optimal_parametersHD
                        )
                        + backgroundHorizontalDetector
                        - peakHorizontalDetector[:, 1]
                    ),
                )
            ),
        )  ## error of the horizontal detector
        fitLevel1_2[f"fitLine_{str(i)}"].create_dataset(
            "errorVerticalDetector",
            dtype="f",
            data=np.transpose(
                (
                    peakVerticalDetector[:, 0],
                    np.absolute(
                        splitPseudoVoigt(
                            peakVerticalDetector[:, 0], optimal_parametersVD
                        )
                        + backgroundVerticalDetector
                        - peakVerticalDetector[:, 1]
                    ),
                )
            ),
        )  ## error of the vertical detector
        for n in range(doublet[i]):
            fitParamsHD = np.append(
                fitParamsHD,
                np.append(
                    optimal_parametersHD[5 * n : 5 * n + 5],
                    [
                        infodictHD["reduced_chisq"],
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
                    ],
                ),
                axis=0,
            )  ##
            fitParamsVD = np.append(
                fitParamsVD,
                np.append(
                    optimal_parametersVD[5 * n : 5 * n + 5],
                    [
                        infodictVD["reduced_chisq"],
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
                    ],
                ),
                axis=0,
            )  ##
            uncertaintyFitParamsHD = np.append(
                uncertaintyFitParamsHD, infodictHD["uncertainties"], axis=0
            )  ##
            uncertaintyFitParamsVD = np.append(
                uncertaintyFitParamsVD, infodictVD["uncertainties"], axis=0
            )  ##

    rawDataLevel1_1.create_dataset(
        "horizontalDetector", dtype="f", data=patternHorizontalDetector
    )  ## save raw data of the horizontal detector
    rawDataLevel1_1.create_dataset(
        "verticalDetector", dtype="f", data=patternVerticalDetector
    )  ## save raw data of the vertical detector
    fitLevel1_2["fitParams"].create_dataset(
        "fitParamsHD",
        dtype="f",
        data=np.reshape(fitParamsHD, (int(np.size(fitParamsHD) / 7), 7)),
    )  ## save parameters of the fit of HD
    fitLevel1_2["fitParams"].create_dataset(
        "fitParamsVD",
        dtype="f",
        data=np.reshape(fitParamsVD, (int(np.size(fitParamsVD) / 7), 7)),
    )  ## save parameters of the fit of VD
    fitLevel1_2["fitParams"].create_dataset(
        "uncertaintyParamsHD", dtype="f", data=uncertaintyFitParamsHD
    )  ## save uncertainty on the parameters of the fit of HD
    fitLevel1_2["fitParams"].create_dataset(
        "uncertaintyParamsVD", dtype="f", data=uncertaintyFitParamsVD
    )  ## save uncertainty on the parameters of the fit of VD

    h5Save.close()
    return
