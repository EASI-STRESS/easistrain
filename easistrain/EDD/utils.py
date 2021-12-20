import argparse
from typing import Callable, Sequence, Tuple, Union
from pathlib import Path
import yaml
import numpy as np
import scipy.optimize
import silx.math.fit


def read_config_file(path: Union[str, Path]):
    with open(path, "r") as config_file:
        return yaml.load(config_file, Loader=yaml.SafeLoader)


def run_from_cli(f: Callable):
    parser = argparse.ArgumentParser()
    parser.add_argument("config_file", type=str, help="Path to the config file")
    args = parser.parse_args()

    config = read_config_file(args.config_file)
    f(**config)
    print("Program finished !")


def linefunc(a, xData):
    return a * xData


def splitPseudoVoigt(xData, *params):
    return silx.math.fit.sum_splitpvoigt(xData, *params)


def gaussEstimation(xData, *params):
    return silx.math.fit.sum_gauss(xData, *params)


def calcBackground(
    xData, yData, fwhmRight, fwhmLeft, counterOfBoxes, nbPeaksInBoxes, guessedPeaksIndex
):
    if int(guessedPeaksIndex[0] - 3 * fwhmLeft) < 0 and int(
        guessedPeaksIndex[-1] + 3 * fwhmRight
    ) <= len(
        xData
    ):  ## case of not enough of points at left
        # print('## case of not enough of points at left')
        xBackground = xData[int(guessedPeaksIndex[-1] + 3 * fwhmRight) :]
        yBackground = yData[int(guessedPeaksIndex[-1] + 3 * fwhmRight) :]
    if (
        int(guessedPeaksIndex[-1] + 3 * fwhmRight) > len(xData)
        and int(guessedPeaksIndex[0] - 3 * fwhmLeft) >= 0
    ):  ## case of not enough of points at right
        # print('## case of not enough of points at right')
        xBackground = xData[0 : int(guessedPeaksIndex[0] - 3 * fwhmLeft)]
        yBackground = yData[0 : int(guessedPeaksIndex[0] - 3 * fwhmLeft)]
    if int(guessedPeaksIndex[0] - 3 * fwhmLeft) < 0 and int(
        guessedPeaksIndex[-1] + 3 * fwhmRight
    ) > len(
        xData
    ):  ## case of not enough of points at left and right
        # print('## case of not enough of points at left and right')
        xBackground = np.append(xData[0:5], xData[-5:])
        yBackground = np.append(yData[0:5], yData[-5:])
    if int(guessedPeaksIndex[0] - 3 * fwhmLeft) >= 0 and int(
        guessedPeaksIndex[-1] + 3 * fwhmRight
    ) <= len(
        xData
    ):  ## case of enough of points at left and right
        # print('## case of enough of points at left and right')
        xBackground = np.append(
            xData[0 : int(guessedPeaksIndex[0] - 3 * fwhmLeft)],
            xData[int(guessedPeaksIndex[-1] + 3 * fwhmRight) :],
        )
        yBackground = np.append(
            yData[0 : int(guessedPeaksIndex[0] - 3 * fwhmLeft)],
            yData[int(guessedPeaksIndex[-1] + 3 * fwhmRight) :],
        )
    # print(xData[0:int(guessedPeaksIndex[0] - 3 * fwhmLeft)])
    # print(xData[-int(guessedPeaksIndex[-1] + 3 * fwhmRight):])
    # print(int(guessedPeaksIndex[-1] + 3 * fwhmRight))
    # print(int(guessedPeaksIndex[0] - 3 * fwhmLeft))
    # print(fwhmRight)
    # print(xBackground)
    # print(yBackground)
    backgroundCoefficient = np.polyfit(
        x=xBackground, y=yBackground, deg=1
    )  ## fit of background with 1d polynom function
    yCalculatedBackground = np.poly1d(backgroundCoefficient)(
        xData
    )  ## yBackground calcuated with the 1d polynom fitted coefficient
    return yCalculatedBackground, backgroundCoefficient


def guessParameters(
    xData: np.ndarray,
    yData: np.ndarray,
    nb_peaks: int,
    withBounds: bool,
):
    fwhmGuess = silx.math.fit.peaks.guess_fwhm(yData)
    peaksGuess = silx.math.fit.peaks.peak_search(
        yData,
        fwhmGuess,
        sensitivity=2.5,
        begin_index=None,
        end_index=None,
        debug=False,
        relevance_info=False,
    )  ## index of the peak with peak relevance

    if np.size(peaksGuess) > nb_peaks:  ## case if more peaks than expected are detected
        peaksGuess = silx.math.fit.peaks.peak_search(
            yData,
            fwhmGuess,
            sensitivity=1,
            begin_index=None,
            end_index=None,
            debug=False,
            relevance_info=True,
        )  ## index of the peak with peak relevance
        peaksGuessArray = np.asarray(peaksGuess)
        orderedIndex = np.argsort(peaksGuessArray[:, 1])[-nb_peaks:]
        peaksGuess = sorted(peaksGuessArray[orderedIndex[:], 0])  ## peaks indices
    if np.size(peaksGuess) < nb_peaks:  ## case if less peaks than expected are detected
        peaksGuess = silx.math.fit.peaks.peak_search(
            yData,
            fwhmGuess,
            sensitivity=1,
            begin_index=None,
            end_index=None,
            debug=False,
            relevance_info=True,
        )  ## index of the peak with peak relevance
        peaksGuessArray = np.asarray(peaksGuess)
        orderedIndex = np.argsort(peaksGuessArray[:, 1])[-nb_peaks:]
        peaksGuess = sorted(peaksGuessArray[orderedIndex[:], 0])  ## peaks indices

    p0Guess = np.zeros(3 * nb_peaks, float)
    minBounds = np.empty((nb_peaks, 3))
    maxBounds = np.empty((nb_peaks, 3))

    for ipar in range(nb_peaks):
        p0Guess[3 * ipar] = yData[int(peaksGuess[ipar])]
        p0Guess[3 * ipar + 1] = xData[int(peaksGuess[ipar])]
        p0Guess[3 * ipar + 2] = fwhmGuess

        minBounds[ipar, 0] = np.amin(yData)  # min H
        minBounds[ipar, 1] = p0Guess[3 * ipar + 1] - 3 * p0Guess[3 * ipar + 2]  # min C
        minBounds[ipar, 2] = 0  # min FWHM

        maxBounds[ipar, 0] = np.amax(yData)  # max H
        maxBounds[ipar, 1] = p0Guess[3 * ipar + 1] + 3 * p0Guess[3 * ipar + 2]  # max C
        maxBounds[ipar, 2] = 2 * p0Guess[3 * ipar + 2]  # max FWHM

    firstGuess, covGuess = scipy.optimize.curve_fit(
        gaussEstimation,
        xData,
        yData,
        p0Guess,
        **(
            {"bounds": (minBounds.flatten(), maxBounds.flatten()), "maxfev": 10000}
            if withBounds
            else {}
        ),
    )

    return firstGuess, peaksGuess


def uChEConversion(a, b, c, ch, ua, ub, uc, uch):
    """
    Calculates the uncertainty on the energy coming from the conversion from channel to energy.

    It includes the uncertainty on the calibration of the coefficients and the peak position
    """
    return np.sqrt(
        ((ua ** 2) * (ch ** 4))
        + ((uch ** 2)) * (((2 * a * ch) + b) ** 2)
        + ((ub ** 2) * ch ** 2)
        + (uc ** 2)
    )


def process_detector_data(fit_min: float, fit_max: float, input_data: np.ndarray):

    peakHorizontalDetector = np.transpose(
        (
            np.arange(fit_min, fit_max),
            input_data[fit_min:fit_max],
        )
    )  ## peak of the horizontal detector
    backgroundHorizontalDetector = silx.math.fit.strip(
        data=peakHorizontalDetector[:, 1],
        w=5,
        niterations=4000,
        factor=1,
        anchors=None,
    )  ## stripped background of the horizontal detector (obtained by stripping the yData)

    return backgroundHorizontalDetector, peakHorizontalDetector
