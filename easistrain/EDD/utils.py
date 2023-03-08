import argparse
from typing import Callable, List, Sequence, Tuple, Union
from pathlib import Path
import yaml
import numpy
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
    print(f"Data processing step '{f.__name__}' finished !")


def linefunc(a, xData):
    return a * xData


def splitPseudoVoigt(xData, *params):
    return silx.math.fit.sum_splitpvoigt(xData, *params)  # type: ignore


def gaussEstimation(xData, *params):
    return silx.math.fit.sum_gauss(xData, *params)  # type: ignore


def calcBackground(
    xData: numpy.ndarray,
    yData: numpy.ndarray,
    fwhmRight: float,
    fwhmLeft: float,
    guessedPeaksIndex: Sequence[float],
) -> numpy.ndarray:
    """Extracts the data outside of the peak bounds to fit the background"""
    try:
        peaks_left_bound = int(guessedPeaksIndex[0] - 3 * fwhmLeft)
        peaks_right_bound = int(guessedPeaksIndex[-1] + 3 * fwhmRight)
        # print(peaks_left_bound, peaks_right_bound)
        if peaks_left_bound <= 5 and peaks_right_bound < len(xData) - 5:
            # case of not enough of points at left
            xBackground = numpy.append(xData[0], xData[peaks_right_bound:])
            yBackground = numpy.append(yData[0], yData[peaks_right_bound:])
        if peaks_right_bound >= len(xData) - 5 and peaks_left_bound > 5:
            # case of not enough of points at right
            xBackground = numpy.append(xData[:peaks_left_bound], xData[-1])
            yBackground = numpy.append(yData[:peaks_left_bound], yData[-1])
        if peaks_left_bound <= 5 and peaks_right_bound >= len(xData) - 5:
            # case of not enough of points at left and right
            xBackground = numpy.append(xData[:5], xData[-5:])
            yBackground = numpy.append(yData[:5], yData[-5:])
        if peaks_left_bound > 5 and peaks_right_bound < len(xData) - 5:
            # case of enough of points at left and right
            xBackground = numpy.append(
                xData[:peaks_left_bound], xData[peaks_right_bound:]
            )
            yBackground = numpy.append(
                yData[:peaks_left_bound], yData[peaks_right_bound:]
            )

        backgroundCoefficient = numpy.polyfit(
            x=xBackground, y=yBackground, deg=1
        )  ## fit of background with 1d polynom function
    except ValueError:
        backgroundCoefficient = numpy.empty(2)
        backgroundCoefficient.fill(numpy.NaN)
    return numpy.poly1d(backgroundCoefficient)(
        xData
    )  ## yBackground calculated with the 1d polynom fitted coefficient


def guessParameters(
    xData: numpy.ndarray,
    yData: numpy.ndarray,
    nb_peaks: int,
    withBounds: bool,
) -> Tuple[numpy.ndarray, List[float]]:
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
        peaksGuessArray = numpy.asarray(raw_peak_guess)
        # Take nb_peaks that are the most relevant
        orderedIndex = numpy.argsort(peaksGuessArray[:, 1])[-nb_peaks:]
        unsorted_peaks_guess: List[float] = peaksGuessArray[
            orderedIndex[:], 0
        ].tolist()  ## peaks indices
    else:
        unsorted_peaks_guess = first_peaks_guess
    peaksGuess = sorted(unsorted_peaks_guess)

    p0Guess = numpy.zeros(3 * nb_peaks, float)
    minBounds = numpy.empty((nb_peaks, 3))
    maxBounds = numpy.empty((nb_peaks, 3))
    try:
        for ipar in range(nb_peaks):
            p0Guess[3 * ipar] = yData[int(peaksGuess[ipar])]
            p0Guess[3 * ipar + 1] = xData[int(peaksGuess[ipar])]
            p0Guess[3 * ipar + 2] = fwhmGuess

            minBounds[ipar, 0] = numpy.amin(yData)  # min H
            minBounds[ipar, 1] = (
                p0Guess[3 * ipar + 1] - 3 * p0Guess[3 * ipar + 2]
            )  # min C
            minBounds[ipar, 2] = 0  # min FWHM

            maxBounds[ipar, 0] = numpy.amax(yData)  # max H
            maxBounds[ipar, 1] = (
                p0Guess[3 * ipar + 1] + 3 * p0Guess[3 * ipar + 2]
            )  # max C
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
    except (IndexError, RuntimeError, ValueError):
        peaksGuess = numpy.empty(nb_peaks)
        peaksGuess.fill(numpy.NaN)
        firstGuess = numpy.empty_like(p0Guess)
        firstGuess.fill(numpy.NaN)
        print("!! guessing fit parameters failed !!")
    return firstGuess, peaksGuess


def uChEConversion(a, b, c, ch, ua, ub, uc, uch):
    """
    Calculates the uncertainty on the energy coming from the conversion from channel to energy.

    It includes the uncertainty on the calibration of the coefficients and the peak position
    """
    return numpy.sqrt(
        ((ua**2) * (ch**4))
        + ((uch**2)) * (((2 * a * ch) + b) ** 2)
        + ((ub**2) * ch**2)
        + (uc**2)
    )


def fit_detector_data(
    channels: numpy.ndarray,
    raw_data: numpy.ndarray,
    nb_peaks: int,
    boxCounter: int,
    scanNumber: int,
    detectorName: str,
):
    """
    Process detector data:
      - Find the background
      - Make a first guess of peak parameters
      - Calculate a background (?)
      - Fit the data without background starting from the first guess
    """
    # print(raw_data, channels)
    first_guess_background: numpy.ndarray = silx.math.fit.strip(  # type: ignore
        data=raw_data,
        w=5,
        niterations=4000,
        factor=1,
        anchors=None,
    )

    peak_guesses, peak_indices = guessParameters(
        xData=channels,
        yData=raw_data - first_guess_background,
        nb_peaks=nb_peaks,
        withBounds=True,
    )  ## guess fit parameters

    calculated_background = calcBackground(
        xData=channels,
        yData=raw_data,
        fwhmRight=peak_guesses[-1],
        fwhmLeft=peak_guesses[2],
        guessedPeaksIndex=peak_indices,
    )

    initial_fit_guess = numpy.zeros(5 * nb_peaks)
    fit_min_bounds = numpy.zeros(5 * nb_peaks)
    fit_max_bounds = numpy.zeros(5 * nb_peaks)
    for n in range(nb_peaks):
        initial_fit_guess[5 * n] = peak_guesses[3 * n]
        initial_fit_guess[5 * n + 1] = peak_guesses[3 * n + 1]
        initial_fit_guess[5 * n + 2] = peak_guesses[3 * n + 2]
        initial_fit_guess[5 * n + 3] = peak_guesses[3 * n + 2]
        initial_fit_guess[5 * n + 4] = 0.5
        fit_min_bounds[5 * n : 5 * n + 5] = [
            numpy.amin(raw_data),
            initial_fit_guess[5 * n + 1] - 3 * initial_fit_guess[5 * n + 2],
            0,
            0,
            0,
        ]
        fit_max_bounds[5 * n : 5 * n + 5] = [
            numpy.amax(raw_data),
            initial_fit_guess[5 * n + 1] + 3 * initial_fit_guess[5 * n + 2],
            len(channels) / 2,
            len(channels) / 2,
            1,
        ]
    try:
        optimal_parameters, covariance = scipy.optimize.curve_fit(
            f=splitPseudoVoigt,
            xdata=channels,
            ydata=raw_data - calculated_background,
            p0=initial_fit_guess,
            bounds=(fit_min_bounds, fit_max_bounds),
            maxfev=10000,
            sigma=numpy.sqrt(0.5 + raw_data),
        )
    except (RuntimeError, ValueError):
        print(
            f"!!Fitting of Peaks in box {boxCounter} in scan {scanNumber} failed for the {detectorName} detector !!"
        )
        print("!! Filling fit parameters with NaN values")
        optimal_parameters = numpy.empty_like(initial_fit_guess)
        optimal_parameters.fill(numpy.NaN)
        covariance = numpy.empty((5 * nb_peaks, 5 * nb_peaks))
        covariance.fill(numpy.NaN)

    fitted_data = splitPseudoVoigt(channels, optimal_parameters) + calculated_background

    goodness_factor = (
        100 * numpy.sum(numpy.absolute(fitted_data - raw_data)) / numpy.sum(raw_data)
    )
    fit_params = numpy.empty((nb_peaks, 6), dtype=optimal_parameters.dtype)
    for n in range(nb_peaks):
        fit_params[n, :5] = optimal_parameters[5 * n : 5 * n + 5]
        fit_params[n, 5] = goodness_factor

    uncertainty_fit_params = numpy.sqrt(numpy.diag(covariance))

    return (
        calculated_background,
        fitted_data,
        fit_params.flatten(),
        uncertainty_fit_params,
    )
