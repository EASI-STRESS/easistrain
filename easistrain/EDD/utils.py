import argparse
from typing import Callable, List, Sequence, Tuple, Union
from pathlib import Path
import yaml
import numpy as np
import scipy.optimize
import silx.math.fit


def read_config_file(path: Union[str, Path]) -> dict:
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
    return silx.math.fit.sum_splitpvoigt(xData, *params)  # type: ignore


def gaussEstimation(xData, *params):
    return silx.math.fit.sum_gauss(xData, *params)  # type: ignore


def calcBackground(
    xData: np.ndarray,
    yData: np.ndarray,
    fwhmRight: float,
    fwhmLeft: float,
    guessedPeaksIndex: Sequence[float],
) -> np.ndarray:
    if int(guessedPeaksIndex[0] - 3 * fwhmLeft) < 0 and int(
        guessedPeaksIndex[-1] + 3 * fwhmRight
    ) <= len(
        xData
    ):  ## case of not enough of points at left
        xBackground = xData[int(guessedPeaksIndex[-1] + 3 * fwhmRight) :]
        yBackground = yData[int(guessedPeaksIndex[-1] + 3 * fwhmRight) :]
    elif (
        int(guessedPeaksIndex[-1] + 3 * fwhmRight) > len(xData)
        and int(guessedPeaksIndex[0] - 3 * fwhmLeft) >= 0
    ):  ## case of not enough of points at right
        xBackground = xData[0 : int(guessedPeaksIndex[0] - 3 * fwhmLeft)]
        yBackground = yData[0 : int(guessedPeaksIndex[0] - 3 * fwhmLeft)]
    elif int(guessedPeaksIndex[0] - 3 * fwhmLeft) < 0 and int(
        guessedPeaksIndex[-1] + 3 * fwhmRight
    ) > len(
        xData
    ):  ## case of not enough of points at left and right
        xBackground = np.append(xData[0:5], xData[-5:])
        yBackground = np.append(yData[0:5], yData[-5:])
    elif int(guessedPeaksIndex[0] - 3 * fwhmLeft) >= 0 and int(
        guessedPeaksIndex[-1] + 3 * fwhmRight
    ) <= len(
        xData
    ):  ## case of enough of points at left and right
        xBackground = np.append(
            xData[0 : int(guessedPeaksIndex[0] - 3 * fwhmLeft)],
            xData[int(guessedPeaksIndex[-1] + 3 * fwhmRight) :],
        )
        yBackground = np.append(
            yData[0 : int(guessedPeaksIndex[0] - 3 * fwhmLeft)],
            yData[int(guessedPeaksIndex[-1] + 3 * fwhmRight) :],
        )
    else:
        raise ValueError(
            "Met a case that should not have been possible when calculating background"
        )

    backgroundCoefficient = np.polyfit(
        x=xBackground, y=yBackground, deg=1
    )  ## fit of background with 1d polynom function
    yCalculatedBackground = np.poly1d(backgroundCoefficient)(
        xData
    )  ## yBackground calcuated with the 1d polynom fitted coefficient
    return yCalculatedBackground


def guessParameters(
    xData: np.ndarray,
    yData: np.ndarray,
    nb_peaks: int,
    withBounds: bool,
) -> Tuple[np.ndarray, List[float]]:
    fwhmGuess = silx.math.fit.peaks.guess_fwhm(yData)
    first_peaks_guess: List[float] = silx.math.fit.peaks.peak_search(
        yData,
        fwhmGuess,
        sensitivity=2.5,
        relevance_info=False,
    )

    if (
        len(first_peaks_guess) != nb_peaks
    ):  ## case if more or less peaks than expected are detected
        raw_peak_guess: List[Tuple[float, float]] = silx.math.fit.peaks.peak_search(
            yData, fwhmGuess, sensitivity=1, relevance_info=True
        )  ## index of the peak with peak relevance
        peaksGuessArray = np.asarray(raw_peak_guess)
        # Take nb_peaks that are the most relevant
        orderedIndex = np.argsort(peaksGuessArray[:, 1])[-nb_peaks:]
        unsorted_peaks_guess: List[float] = peaksGuessArray[
            orderedIndex[:], 0
        ].tolist()  ## peaks indices
    else:
        unsorted_peaks_guess = first_peaks_guess
    peaksGuess = sorted(unsorted_peaks_guess)

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


def process_detector_data(
    fit_min: float, fit_max: float, input_data: np.ndarray, nb_peaks: int
):
    """
    Process detector data:
      - Find the background
      - Make a first guess of peak parameters
      - Calculate a background (?)
      - Fit the data without background starting from the first guess
    """

    channels: np.ndarray[int, np.float64] = np.arange(fit_min, fit_max)
    raw_data: np.ndarray = input_data[fit_min:fit_max]

    data_background: np.ndarray = silx.math.fit.strip(  # type: ignore
        data=raw_data,
        w=5,
        niterations=4000,
        factor=1,
        anchors=None,
    )  ## background of the horizontal detector (obtained by stripping the yData)

    peak_guesses, peak_indices = guessParameters(
        xData=channels,
        yData=raw_data - data_background,
        nb_peaks=nb_peaks,
        withBounds=True,
    )  ## guess fit parameters for HD

    calculated_background = calcBackground(
        xData=channels,
        yData=raw_data,
        fwhmRight=peak_guesses[-1],
        fwhmLeft=peak_guesses[2],
        guessedPeaksIndex=peak_indices,
    )

    initial_fit_guess = np.zeros(5 * nb_peaks)
    fit_min_bounds = np.zeros(5 * nb_peaks)
    fit_max_bounds = np.zeros(5 * nb_peaks)
    for n in range(nb_peaks):
        initial_fit_guess[5 * n] = peak_guesses[3 * n]
        initial_fit_guess[5 * n + 1] = peak_guesses[3 * n + 1]
        initial_fit_guess[5 * n + 2] = peak_guesses[3 * n + 2]
        initial_fit_guess[5 * n + 3] = peak_guesses[3 * n + 2]
        initial_fit_guess[5 * n + 4] = 0.5
        fit_min_bounds[5 * n : 5 * n + 5] = [
            np.amin(raw_data),
            initial_fit_guess[5 * n + 1] - 3 * initial_fit_guess[5 * n + 2],
            0,
            0,
            0,
        ]  # minimum bounds of the parametrs solution (H, C, FWHM1, FWHM2, eta) to append for the horizontal detector
        fit_max_bounds[5 * n : 5 * n + 5] = [
            np.amax(raw_data),
            initial_fit_guess[5 * n + 1] + 3 * initial_fit_guess[5 * n + 2],
            len(channels) / 2,
            len(channels) / 2,
            1,
        ]
    optimal_parameters, covariance = scipy.optimize.curve_fit(
        f=splitPseudoVoigt,
        xdata=channels,
        ydata=raw_data - calculated_background,
        p0=initial_fit_guess,
        bounds=(fit_min_bounds, fit_max_bounds),
        maxfev=10000,
    )  ## fit of the peak of the Horizontal detector

    fit_params = np.array(())
    uncertainty_fit_params = np.sqrt(np.diag(covariance))
    for n in range(nb_peaks):
        fit_params = np.append(
            fit_params,
            np.append(
                optimal_parameters[5 * n : 5 * n + 5],
                100
                * np.sum(
                    np.absolute(
                        splitPseudoVoigt(
                            channels,
                            optimal_parameters,
                        )
                        + data_background
                        - raw_data
                    )
                )
                / np.sum(raw_data),
            ),
            axis=0,
        )

    return (
        np.transpose(
            (
                channels,
                raw_data,
            )
        ),  # peakHorizontalDetector
        calculated_background,
        optimal_parameters,
        fit_params,
        uncertainty_fit_params,
    )
