import os
from typing import Sequence

import h5py
import numpy

from easistrain.EDD.utils import run_from_cli


def transformationMatrix(
    rx: float, ry: float, rz: float, tx: float, ty: float, tz: float
):
    rx = numpy.deg2rad(rx)
    ry = numpy.deg2rad(ry)
    rz = numpy.deg2rad(rz)
    transfMat = numpy.zeros((4, 4), "float64")
    rotx = numpy.array(
        [
            [1, 0, 0],
            [0, numpy.cos(rx), -numpy.sin(rx)],
            [0, numpy.sin(rx), numpy.cos(rx)],
        ]
    )  ### Rotation matrix around x, x is in the direction of the beam ###
    roty = numpy.array(
        [
            [numpy.cos(ry), 0, numpy.sin(ry)],
            [0, 1, 0],
            [-numpy.sin(ry), 0, numpy.cos(ry)],
        ]
    )  ### Rotation matrix around y, y is perpendicular to the beam and in its plane ###
    rotz = numpy.array(
        [
            [numpy.cos(rz), -numpy.sin(rz), 0],
            [numpy.sin(rz), numpy.cos(rz), 0],
            [0, 0, 1],
        ]
    )  ### Rotation matrix around z, z iz perpendicular to the beam and out of the beam plane ###
    rotxyz = numpy.dot(
        rotz, numpy.dot(roty, rotx)
    )  ### Rotation matrix for 3 rotations: first around x, second around y and third around z  ###
    transfMat[0:3, 0:3] = rotxyz[0:3, 0:3]
    transfMat[3, 0:3] = 0
    transfMat[0, 3] = tx
    transfMat[1, 3] = ty
    transfMat[2, 3] = tz
    transfMat[3, 3] = 1
    return transfMat


def coordTransformation(
    fileRead: str, fileSave: str, numberOfPeaks: int, gonioToSample: Sequence[float]
):
    with h5py.File(fileRead, "r") as h5Read:  ## Read the h5 file of raw data
        scanList = list(h5Read.keys())  ## list of the scans
        lengthCounter = 0
        for scan in scanList:
            shapeDsetPeak = h5Read[f"{scan}/tthPositionsGroup/peak_0000"].shape[0]
            lengthCounter = (
                lengthCounter + shapeDsetPeak
            )  ## shape of the matrix on which all the points for one peak will be saved
    transfMat = transformationMatrix(
        gonioToSample[0],
        gonioToSample[1],
        gonioToSample[2],
        gonioToSample[3],
        gonioToSample[4],
        gonioToSample[5],
    )

    if os.path.dirname(fileSave):
        os.makedirs(os.path.dirname(fileSave), exist_ok=True)
    with h5py.File(fileSave, "a") as h5Save:  ## create/append h5 file to save in
        globalGroup = h5Save.create_group(
            "global",
        )  ## Creation of the global group in which all peaks positions of all the points will be put

        for peakNumber in range(numberOfPeaks):
            globalPeak = globalGroup.create_dataset(
                f"peak_{str(peakNumber).zfill(4)}",
                dtype="float64",
                data=numpy.zeros((lengthCounter, 13), "float64"),
            )  ## creation of the dataset of peak info (coordinate in gonio, center, ...) for each peak
            globalPeakInSample = globalGroup.create_dataset(
                f"inSample_peak_{str(peakNumber).zfill(4)}",
                dtype="float64",
                data=numpy.zeros((lengthCounter, 13), "float64"),
            )  ## creation of the dataset of peak info (coordinate in sample, center, ...) for each peak
            uncertaintyGlobalPeak = globalGroup.create_dataset(
                f"uncertaintyPeak_{str(peakNumber).zfill(4)}",
                dtype="float64",
                data=numpy.zeros((lengthCounter, 13), "float64"),
            )  ## creation of the dataset of uncertainty for each peak (gonio coordinates)
            uncertaintyGlobalPeakInSample = globalGroup.create_dataset(
                f"inSample_uncertaintyPeak_{str(peakNumber).zfill(4)}",
                dtype="float64",
                data=numpy.zeros((lengthCounter, 13), "float64"),
            )  ## creation of the dataset of uncertainty for each peak (sample coordinates)
            rowsCounter = 0
            for scan in scanList:
                with h5py.File(
                    fileRead, "r"
                ) as h5Read:  ## Read the h5 file of raw data
                    shapeDsetPeak = h5Read[f"{scan}/tthPositionsGroup/peak_0000"].shape[
                        0
                    ]
                    globalPeak[rowsCounter : rowsCounter + shapeDsetPeak, :] = h5Read[
                        f"{scan}/tthPositionsGroup/peak_{str(peakNumber).zfill(4)}"
                    ][
                        ()
                    ]  ## put peak info (coordinate in gonio, center, ...) in global dataset regrouping all point for one peak
                    uncertaintyGlobalPeak[
                        rowsCounter : rowsCounter + shapeDsetPeak, :
                    ] = h5Read[
                        f"{scan}/tthPositionsGroup/uncertaintyPeak_{str(peakNumber).zfill(4)}"
                    ][
                        ()
                    ]  ## put uncertainty of peak info (coordinate, center, ...) in global dataset regrouping all point for one peak
                    globalPeakInSample[
                        rowsCounter : rowsCounter + shapeDsetPeak, :
                    ] = h5Read[
                        f"{scan}/tthPositionsGroup/peak_{str(peakNumber).zfill(4)}"
                    ][
                        ()
                    ]  ## same as above, then belwo convert the coordinate from gonio to sample
                    uncertaintyGlobalPeakInSample[
                        rowsCounter : rowsCounter + shapeDsetPeak, :
                    ] = h5Read[
                        f"{scan}/tthPositionsGroup/uncertaintyPeak_{str(peakNumber).zfill(4)}"
                    ][
                        ()
                    ]  ## same as above, then belwo convert the coordinate from gonio to sample
                    rowsCounter = rowsCounter + shapeDsetPeak
            for icount in range(len(globalPeakInSample)):
                coordInSample = numpy.dot(
                    transfMat,
                    [
                        -globalPeakInSample[icount, 0],
                        -globalPeakInSample[icount, 1],
                        -globalPeakInSample[icount, 2],
                        1,
                    ],
                )
                globalPeakInSample[icount, 0] = coordInSample[0]
                globalPeakInSample[icount, 1] = coordInSample[1]
                globalPeakInSample[icount, 2] = coordInSample[2]
            uncertaintyGlobalPeakInSample[:, 0:3] = globalPeakInSample[:, 0:3]


if __name__ == "__main__":
    run_from_cli(coordTransformation)
