import numpy as np
import h5py
import matplotlib.pyplot as plt


def plot1D(
    fileRead,
    sample,
    dataset,
    scanNumber,
    counterBeforeSample,
    counterAfterSample,
    positioner,
):
    xData = np.array(())
    yData = np.array(())
    secXData = np.array(())
    with h5py.File(
        fileRead, "r"
    ) as h5Read:  ## Read the h5 file of the energy calibration of the detectors
        for scanNumber in range(scanNumber[0], scanNumber[1] + 1):
            xDataUpdate = h5Read[
                f"{sample}_{str(dataset)}_{str(scanNumber)}.1/instrument/positioners/{positioner}"
            ][
                ()
            ]  ## the position of the motor to put in the x axis
            yDataUpdate = (
                h5Read[
                    f"{sample}_{str(dataset)}_{str(scanNumber)}.1/measurement/{counterAfterSample}"
                ][()]
                / h5Read[
                    f"{sample}_{str(dataset)}_{str(scanNumber)}.1/measurement/{counterBeforeSample}"
                ][()]
            )  ## The data in the counter to put in the y axis
            xData = np.append(xData, xDataUpdate)
            yData = np.append(yData, yDataUpdate)
            secXData = np.append(
                secXData, scanNumber
            )  ## The scanNumber to put in the axis of the second figure
    plt.figure(figsize=(10, 8))
    plt.plot(xData, yData, ".")
    plt.xlabel(positioner, family="sans-serif", fontsize=28)
    plt.ylabel(
        f"{counterAfterSample}/{counterBeforeSample}", family="sans-serif", fontsize=28
    )
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid()
    plt.figure(figsize=(10, 8))
    plt.plot(secXData, yData, ".")
    plt.xlabel("Scan Number", family="sans-serif", fontsize=28)
    plt.ylabel(
        f"{counterAfterSample}/{counterBeforeSample}", family="sans-serif", fontsize=28
    )
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid()
    plt.show()
    return
