import os
from typing import Sequence

import numpy
import h5py

from easistrain.EDD.constants import pCstInkeVS, speedLightInAPerS
from easistrain.EDD.math import compute_qs
from easistrain.EDD.utils import run_from_cli, uChEConversion


def uE(dspacing, ttheta, udspacing, uttheta):
    """Calculates the uncertainty on the energy coming from d and the fixed diffraction angle"""
    theta = numpy.radians(0.5 * ttheta)
    csctheta = 1 / numpy.sin(theta)
    return numpy.sqrt(
        ((38.4276 * (udspacing**2) * (csctheta**2)) / (dspacing**4))
        + (
            (38.4276 * (uttheta**2) * ((1 / numpy.tan(theta)) ** 2) * (csctheta**2))
            / (dspacing**2)
        )
    )


def ud(energy, ttheta, uenergy, uttheta):
    """Calculates the uncertainty on the dspacing coming from the energy and the fixed diffraction angle"""
    theta = numpy.radians(0.5 * ttheta)
    csctheta = 1 / numpy.sin(theta)
    return numpy.sqrt(
        ((38.4276 * (uenergy**2) * (csctheta**2)) / (energy**4))
        + (
            (38.4276 * (uttheta**2) * ((1 / numpy.tan(theta)) ** 2) * (csctheta**2))
            / (energy**2)
        )
    )


def ustrain(energy0, energystrained, uenergy0, uenergystrained):
    """Calculates the uncertainty on the strain coming from the measured strain-free energy and the strained energy"""
    return numpy.sqrt(
        ((uenergy0**2) / (energy0**2))
        + ((uenergystrained**2) / (energystrained**2))
    )


def _energy_wavelength(x):
    """keV to anstrom and vice versa"""
    return pCstInkeVS * speedLightInAPerS / x


def _channels_to_dspacing(channels, energy_calib_coeff, two_theta):
    energies = numpy.polyval(energy_calib_coeff, channels)
    wavelength = _energy_wavelength(energies)
    dspacing = wavelength / (2 * numpy.sin(numpy.deg2rad(two_theta / 2)))
    return dspacing


def preStraind0cstEDD(
    fileRead: str,
    fileSave: str,
    pathFileDetectorCalibration: str,
    scanDetectorCalibration: str,
    pathFileAngleCalibration: str,
    scanAngleCalibration: str,
    pathFileReferenceFitEDD: str,
    scanReferenceFitEDD: str,
    numberOfPeaks: int,
    d0: Sequence[float],
):
    if os.path.dirname(fileSave):
        os.makedirs(os.path.dirname(fileSave), exist_ok=True)
    with h5py.File(fileSave, "a") as h5Save:  ## create/append h5 file to save in
        strainGroupWithd0 = h5Save.create_group(
            "STRAIN_with_d0"
        )  ## creation of the strain group for results with d0
        # strainGroupWithoutd0 = h5Save.create_group('STRAIN_without_d0') ## creation of the strain group for results without d0

        with h5py.File(
            pathFileDetectorCalibration, "r"
        ) as h5CalibRead:  ## Read the h5 file of the energy calibration of the detectors
            calibCoeffsHD = h5CalibRead[
                f"detectorCalibration/{scanDetectorCalibration}/calibCoeffs/calibCoeffsHD"
            ][
                ()
            ]  ## import the energy calibration coefficients of the horizontal detector
            calibCoeffsVD = h5CalibRead[
                f"detectorCalibration/{scanDetectorCalibration}/calibCoeffs/calibCoeffsVD"
            ][
                ()
            ]  ## import the energy calibration coefficients of the vertical detector
            uncertaintyCalibCoeffsHD = h5CalibRead[
                f"detectorCalibration/{scanDetectorCalibration}/calibCoeffs/uncertaintyCalibCoeffsHD"
            ][
                ()
            ]  ## import the uncertainty of the energy calibration coefficients of the horizontal detector
            uncertaintyCalibCoeffsVD = h5CalibRead[
                f"detectorCalibration/{scanDetectorCalibration}/calibCoeffs/uncertaintyCalibCoeffsVD"
            ][
                ()
            ]  ## import the uncertainty of the energy calibration coefficients of the vertical detector

        with h5py.File(
            pathFileAngleCalibration, "r"
        ) as h5AngleRead:  ## Read the h5 file of the angle calibration of the detectors
            AngleHD = h5AngleRead[
                f"angleCalibration/{scanAngleCalibration}/calibratedAngle/calibratedAngleHD"
            ][
                ()
            ]  ## import the calibrated angle of the horizontal detector
            AngleVD = h5AngleRead[
                f"angleCalibration/{scanAngleCalibration}/calibratedAngle/calibratedAngleVD"
            ][
                ()
            ]  ## import the calibrated angle of the vertical detector

        if pathFileReferenceFitEDD:
            with h5py.File(pathFileReferenceFitEDD, "r") as h5RefFitRead:
                RefChanHD = h5RefFitRead[
                    f"/{scanReferenceFitEDD}/fit/0000/fitParams/fitParamsHD"
                ][:, 1]
                RefChanVD = h5RefFitRead[
                    f"/{scanReferenceFitEDD}/fit/0000/fitParams/fitParamsVD"
                ][:, 1]
                d0HD = _channels_to_dspacing(RefChanHD, calibCoeffsHD, AngleHD)
                d0VD = _channels_to_dspacing(RefChanVD, calibCoeffsVD, AngleVD)
                d0 = (d0HD + d0VD) / 2
                numberOfPeaks = len(d0)

        with h5py.File(fileRead, "r") as h5Read:  ## Read the h5 file of raw data
            for peakNumber in range(numberOfPeaks):
                strainPerPeakWithd0 = strainGroupWithd0.create_group(
                    f"peak_{str(peakNumber).zfill(4)}"
                )  ## create group for each peak in the strain group
                for i in range(
                    int(
                        len(list(h5Read[f"pointsPerPeak_{str(peakNumber).zfill(4)}"]))
                        / 2
                    )
                ):
                    allPtsInPeak = h5Read[
                        f"pointsPerPeak_{str(peakNumber).zfill(4)}/point_{str(i).zfill(5)}"
                    ][
                        ()
                    ]  ##
                    shapeallPtsInPeak = h5Read[
                        f"pointsPerPeak_{str(peakNumber).zfill(4)}/point_{str(i).zfill(5)}"
                    ].shape[
                        0
                    ]  ##
                    uallPtsInPeak = h5Read[
                        f"pointsPerPeak_{str(peakNumber).zfill(4)}/uncertaintyPoint_{str(i).zfill(5)}"
                    ][
                        ()
                    ]  ##
                    pts = strainPerPeakWithd0.create_dataset(
                        f"point_{str(i).zfill(5)}",
                        dtype="float64",
                        data=numpy.zeros((shapeallPtsInPeak, 12), "float64"),
                    )  ## dataset for each point
                    uncertaintyPts = strainPerPeakWithd0.create_dataset(
                        f"uncertainty_point_{str(i).zfill(5)}",
                        dtype="float64",
                        data=numpy.zeros((shapeallPtsInPeak, 12), "float64"),
                    )  ## dataset for uncertainty for each point
                    pts[:, 0:8] = allPtsInPeak[
                        :, 0:8
                    ]  ## Coordinates of the point, goniometrtic angles and beam direction
                    uncertaintyPts[:, 0:8] = uallPtsInPeak[
                        :, 0:8
                    ]  ## Coordinates of the point, goniometrtic angles and beam direction
                    for j in range(shapeallPtsInPeak):
                        if (
                            allPtsInPeak[j, 6] == -90
                            and numpy.polyval(calibCoeffsHD, allPtsInPeak[j, 8]) > 0
                        ):  ## it is the horizontal detector
                            pts[j, 8] = numpy.log(
                                (
                                    (pCstInkeVS * speedLightInAPerS)
                                    / (
                                        2
                                        * d0[peakNumber]
                                        * numpy.sin(numpy.deg2rad(0.5 * AngleHD))
                                    )
                                )
                                / numpy.polyval(calibCoeffsHD, allPtsInPeak[j, 8])
                            )  ## strain measured using the HD
                            uncertaintyPts[j, 8] = ustrain(
                                (
                                    (pCstInkeVS * speedLightInAPerS)
                                    / (
                                        2
                                        * d0[peakNumber]
                                        * numpy.sin(numpy.deg2rad(0.5 * AngleHD))
                                    )
                                ),
                                numpy.polyval(calibCoeffsHD, allPtsInPeak[j, 8]),
                                uChEConversion(
                                    calibCoeffsHD[0],
                                    calibCoeffsHD[1],
                                    calibCoeffsHD[2],
                                    allPtsInPeak[j, 8],
                                    uncertaintyCalibCoeffsHD[0],
                                    uncertaintyCalibCoeffsHD[1],
                                    uncertaintyCalibCoeffsHD[2],
                                    uallPtsInPeak[j, 8],
                                ),
                                uChEConversion(
                                    calibCoeffsHD[0],
                                    calibCoeffsHD[1],
                                    calibCoeffsHD[2],
                                    allPtsInPeak[j, 8],
                                    uncertaintyCalibCoeffsHD[0],
                                    uncertaintyCalibCoeffsHD[1],
                                    uncertaintyCalibCoeffsHD[2],
                                    uallPtsInPeak[j, 8],
                                ),
                            )  ## uncertainty on the strain of the HD
                        if (
                            allPtsInPeak[j, 6] == 0
                            and numpy.polyval(calibCoeffsVD, allPtsInPeak[j, 8]) > 0
                        ):  ## it is the vertical detector
                            pts[j, 8] = numpy.log(
                                (
                                    (pCstInkeVS * speedLightInAPerS)
                                    / (
                                        2
                                        * d0[peakNumber]
                                        * numpy.sin(numpy.deg2rad(0.5 * AngleVD))
                                    )
                                )
                                / numpy.polyval(calibCoeffsVD, allPtsInPeak[j, 8])
                            )  ## strain measured using the VD
                            uncertaintyPts[j, 8] = ustrain(
                                (
                                    (pCstInkeVS * speedLightInAPerS)
                                    / (
                                        2
                                        * d0[peakNumber]
                                        * numpy.sin(numpy.deg2rad(0.5 * AngleVD))
                                    )
                                ),
                                numpy.polyval(calibCoeffsVD, allPtsInPeak[j, 8]),
                                uChEConversion(
                                    calibCoeffsVD[0],
                                    calibCoeffsVD[1],
                                    calibCoeffsVD[2],
                                    allPtsInPeak[j, 8],
                                    uncertaintyCalibCoeffsVD[0],
                                    uncertaintyCalibCoeffsVD[1],
                                    uncertaintyCalibCoeffsVD[2],
                                    uallPtsInPeak[j, 8],
                                ),
                                uChEConversion(
                                    calibCoeffsVD[0],
                                    calibCoeffsVD[1],
                                    calibCoeffsVD[2],
                                    allPtsInPeak[j, 8],
                                    uncertaintyCalibCoeffsVD[0],
                                    uncertaintyCalibCoeffsVD[1],
                                    uncertaintyCalibCoeffsVD[2],
                                    uallPtsInPeak[j, 8],
                                ),
                            )  ## uncertainty on the strain of the VD
                        # print(allPtsInPeak[j, 3:8])
                        qq1, qq2, qq3 = compute_qs(allPtsInPeak[j, 3:8])
                        pts[
                            j, 9
                        ] = qq1  ## component of the scattering vector in the x direction
                        pts[
                            j, 10
                        ] = qq2  ## component of the scattering vector in the y direction
                        pts[
                            j, 11
                        ] = qq3  ## component of the scattering vector in the z direction
                        uncertaintyPts[
                            j, 9
                        ] = qq1  ## component of the scattering vector in the x direction
                        uncertaintyPts[
                            j, 10
                        ] = qq2  ## component of the scattering vector in the y direction
                        uncertaintyPts[
                            j, 11
                        ] = qq3  ## component of the scattering vector in the z direction


if __name__ == "__main__":
    run_from_cli(preStraind0cstEDD)
