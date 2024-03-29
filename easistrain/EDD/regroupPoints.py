import os
from typing import Sequence

import numpy
import h5py

from easistrain.EDD.utils import run_from_cli


def regroupPoints(fileRead: Sequence[str], fileSave: str, numberOfPeaks: int):
    if os.path.dirname(fileSave):
        os.makedirs(os.path.dirname(fileSave), exist_ok=True)
    with h5py.File(fileSave, "a") as h5Save:  ## create/append h5 file to save in
        rowsInAll = 0
        for fileR in fileRead:
            with h5py.File(fileR, "r") as h5Read:  ## Read the h5 file of raw data
                countRows = h5Read["global/inSample_peak_0000"].shape[0]  ## shape of
                # prvrowsInAll = rowsInAll
                rowsInAll = rowsInAll + countRows
        countFiller = 0
        for fileR in fileRead:
            for peakNumber in range(numberOfPeaks):
                if (
                    f"coordInSample_Peak_{str(peakNumber).zfill(4)}"
                    not in h5Save.keys()
                ):
                    pointsInPeakGlobal = h5Save.create_dataset(
                        f"coordInSample_Peak_{str(peakNumber).zfill(4)}",
                        dtype="float64",
                        data=numpy.zeros((rowsInAll, 13), "float64"),
                    )
                else:
                    pointsInPeakGlobal = h5Save[
                        f"coordInSample_Peak_{str(peakNumber).zfill(4)}"
                    ]
                if (
                    f"coordInSample_uncertainty_Peak_{str(peakNumber).zfill(4)}"
                    not in h5Save.keys()
                ):
                    upointsInPeakGlobal = h5Save.create_dataset(
                        f"coordInSample_uncertainty_Peak_{str(peakNumber).zfill(4)}",
                        dtype="float64",
                        data=numpy.zeros((rowsInAll, 13), "float64"),
                    )
                else:
                    upointsInPeakGlobal = h5Save[
                        f"coordInSample_uncertainty_Peak_{str(peakNumber).zfill(4)}"
                    ]
                with h5py.File(fileR, "r") as h5Read:  ## Read the h5 file of raw data
                    shapeScanPoints = h5Read[
                        f"global/inSample_peak_{str(peakNumber).zfill(4)}"
                    ].shape[
                        0
                    ]  ## length
                    pointsInPeakGlobal[
                        countFiller : countFiller + shapeScanPoints, :
                    ] = h5Read[f"global/inSample_peak_{str(peakNumber).zfill(4)}"][
                        ()
                    ]  ## points of the peak
                    upointsInPeakGlobal[
                        countFiller : countFiller + shapeScanPoints, :
                    ] = h5Read[
                        f"global/inSample_uncertaintyPeak_{str(peakNumber).zfill(4)}"
                    ][
                        ()
                    ]  ## uncertainty of the points of the peak
            countFiller = countFiller + shapeScanPoints
        for peakNumber in range(numberOfPeaks):
            peakGroup = h5Save.create_group(
                f"pointsPerPeak_{str(peakNumber).zfill(4)}"
            )  ## create a group per peak in the global group
            pointsCounter = 0
            matToFilter = h5Save[f"coordInSample_Peak_{str(peakNumber).zfill(4)}"][
                ()
            ]  ## the peak matrix to filter
            umatToFilter = h5Save[
                f"coordInSample_uncertainty_Peak_{str(peakNumber).zfill(4)}"
            ][
                ()
            ]  ## the uncertainty matrix to filter
            matCheck = numpy.array([[]])
            umatCheck = numpy.array([[]])
            for i in range(len(matToFilter)):
                matFirstFilter = matToFilter[
                    numpy.round(matToFilter[:, 0], 4)
                    == numpy.round(matToFilter[i, 0], 4),
                    :,
                ]  ## filtering of peak info based on x coordinate of the measurement point
                matSecFilter = matFirstFilter[
                    numpy.round(matFirstFilter[:, 1], 4)
                    == numpy.round(matToFilter[i, 1], 4),
                    :,
                ]  ## filtering of peak info based on y coordinate of the measurement point
                matThirdFilter = matSecFilter[
                    numpy.round(matSecFilter[:, 2], 4)
                    == numpy.round(matToFilter[i, 2], 4),
                    :,
                ]  ## filtering of peak info based on z coordinate of the measurement point
                umatFirstFilter = umatToFilter[
                    numpy.round(umatToFilter[:, 0], 4)
                    == numpy.round(umatToFilter[i, 0], 4),
                    :,
                ]  ## filtering of uncertainty based on x coordinate of the measurement point
                umatSecFilter = umatFirstFilter[
                    numpy.round(umatFirstFilter[:, 1], 4)
                    == numpy.round(umatToFilter[i, 1], 4),
                    :,
                ]  ## filtering of uncertainty based on y coordinate of the measurement point
                umatThirdFilter = umatSecFilter[
                    numpy.round(umatSecFilter[:, 2], 4)
                    == numpy.round(umatToFilter[i, 2], 4),
                    :,
                ]  ## filtering of uncertainty based on z coordinate of the measurement point
                # print(matThirdFilter)
                if pointsCounter == 0:
                    peakGroup.create_dataset(
                        f"point_{str(pointsCounter).zfill(5)}",
                        dtype="float64",
                        data=matThirdFilter,
                    )  ## save of peak info after separation based on point coordinates (! just for the first scan)
                    peakGroup.create_dataset(
                        f"uncertaintyPoint_{str(pointsCounter).zfill(5)}",
                        dtype="float64",
                        data=umatThirdFilter,
                    )  ## save of uncertainty after separation based on point coordinates (! just for the first scan)
                    matCheck = numpy.append(matCheck, matThirdFilter[0, :3])
                    matCheck = numpy.reshape(matCheck, (int(len(matCheck) / 3), 3))
                    umatCheck = numpy.append(umatCheck, umatThirdFilter[0, :3])
                    umatCheck = numpy.reshape(umatCheck, (int(len(umatCheck) / 3), 3))
                    pointsCounter = pointsCounter + 1
                if pointsCounter > 0:
                    check = numpy.array(())
                    for kk in range(len(matCheck)):
                        check = numpy.append(
                            check, numpy.sum(matThirdFilter[0, :3] == matCheck[kk, :3])
                        )
                    if numpy.sum(check[:] == 3) == 0:
                        peakGroup.create_dataset(
                            f"point_{str(pointsCounter).zfill(5)}",
                            dtype="float64",
                            data=matThirdFilter,
                        )  ## save of peak info after separation based on point coordinates (! for the rest of the scans)
                        peakGroup.create_dataset(
                            f"uncertaintyPoint_{str(pointsCounter).zfill(5)}",
                            dtype="float64",
                            data=umatThirdFilter,
                        )  ## save of uncertainty after separation based on point coordinates (! for the rest of the scans)
                        pointsCounter = pointsCounter + 1
                        matCheck = numpy.append(matCheck, matThirdFilter[0, :3])
                        matCheck = numpy.reshape(matCheck, (int(len(matCheck) / 3), 3))
                        umatCheck = numpy.append(umatCheck, umatThirdFilter[0, :3])
                        umatCheck = numpy.reshape(
                            umatCheck, (int(len(umatCheck) / 3), 3)
                        )
                    else:
                        pointsCounter = pointsCounter


if __name__ == "__main__":
    run_from_cli(regroupPoints)
