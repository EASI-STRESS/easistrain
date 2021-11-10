import argparse
from typing import Callable, Union
from pathlib import Path
import yaml
import numpy as np
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
