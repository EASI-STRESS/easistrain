import numpy as np
import h5py

from easistrain.EDD.constants import pCstInkeVS, speedLightInAPerS
from easistrain.EDD.utils import run_from_cli, uChEConversion


def uE(dspacing, ttheta, udspacing, uttheta):
    """Calculates the uncertainty on the energy coming from d and the fixed diffraction angle"""
    theta = np.radians(0.5 * ttheta)
    csctheta = 1 / np.sin(theta)
    return np.sqrt(
        ((38.4276 * (udspacing ** 2) * (csctheta ** 2)) / (dspacing ** 4))
        + (
            (38.4276 * (uttheta ** 2) * ((1 / np.tan(theta)) ** 2) * (csctheta ** 2))
            / (dspacing ** 2)
        )
    )


def ud(energy, ttheta, uenergy, uttheta):
    """Calculates the uncertainty on the dspacing coming from the energy and the fixed diffraction angle"""
    theta = np.radians(0.5 * ttheta)
    csctheta = 1 / np.sin(theta)
    return np.sqrt(
        ((38.4276 * (uenergy ** 2) * (csctheta ** 2)) / (energy ** 4))
        + (
            (38.4276 * (uttheta ** 2) * ((1 / np.tan(theta)) ** 2) * (csctheta ** 2))
            / (energy ** 2)
        )
    )


def ustrain(energy0, energystrained, uenergy0, uenergystrained):
    """Calculates the uncertainty on the strain coming from the measured strain-free energy and the strained energy"""
    return np.sqrt(
        ((uenergy0 ** 2) / (energy0 ** 2))
        + ((uenergystrained ** 2) / (energystrained ** 2))
    )


def diffVector(angles):
    """Calcuates the direction of the q vector in the sample reference space using the goniometric angles, diffraction angle and azimuth angle"""
    phi = np.radians(angles[0])
    chi = np.radians(angles[1])
    omega = np.radians(angles[2])
    delta = np.radians(angles[3])
    theta = np.radians(0.5 * angles[4])
    q1 = (
        (np.cos(theta) * np.cos(chi) * np.sin(delta) * np.sin(phi))
        + (
            np.cos(delta)
            * np.cos(theta)
            * (
                (np.cos(phi) * np.sin(omega))
                - (np.cos(omega) * np.sin(phi) * np.sin(chi))
            )
        )
        - np.sin(theta)
        * ((np.cos(phi) * np.cos(omega)) + (np.sin(phi) * np.sin(chi) * np.sin(omega)))
    )
    q2 = (
        np.cos(delta)
        * np.cos(theta)
        * ((np.cos(phi) * np.cos(omega) * np.sin(chi)) + (np.sin(phi) * np.sin(omega)))
        - (np.cos(theta) * np.cos(phi) * np.cos(chi) * np.sin(delta))
        - (
            np.sin(theta)
            * (
                (np.cos(omega) * np.sin(phi))
                - (np.cos(phi) * np.sin(chi) * np.sin(omega))
            )
        )
    )
    q3 = (
        (np.cos(delta) * np.cos(theta) * np.cos(chi) * np.cos(omega))
        + (np.cos(theta) * np.sin(delta) * np.sin(chi))
        + (np.cos(chi) * np.sin(theta) * np.sin(omega))
    )
    return (
        q1,
        q2,
        q3,
    )


def preStraind0cstEDD(
    fileRead,
    fileSave,
    pathFileDetectorCalibration,
    scanDetectorCalibration,
    pathFileAngleCalibration,
    scanAngleCalibration,
    numberOfPeaks,
    d0,
):

    h5Save = h5py.File(fileSave, "a")  ## create/append h5 file to save in
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

    with h5py.File(fileRead, "r") as h5Read:  ## Read the h5 file of raw data
        for peakNumber in range(numberOfPeaks):
            strainPerPeakWithd0 = strainGroupWithd0.create_group(
                f"peak_{str(peakNumber).zfill(4)}"
            )  ## create group for each peak in the strain group
            for i in range(
                int(len(list(h5Read[f"pointsPerPeak_{str(peakNumber).zfill(4)}"])) / 2)
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
                    data=np.zeros((shapeallPtsInPeak, 12), "float64"),
                )  ## dataset for each point
                uncertaintyPts = strainPerPeakWithd0.create_dataset(
                    f"uncertainty_point_{str(i).zfill(5)}",
                    dtype="float64",
                    data=np.zeros((shapeallPtsInPeak, 12), "float64"),
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
                        and np.polyval(calibCoeffsHD, allPtsInPeak[j, 8]) > 0
                    ):  ## it is the horizontal detector
                        pts[j, 8] = np.log(
                            (
                                (pCstInkeVS * speedLightInAPerS)
                                / (
                                    2
                                    * d0[peakNumber]
                                    * np.sin(np.deg2rad(0.5 * AngleHD))
                                )
                            )
                            / np.polyval(calibCoeffsHD, allPtsInPeak[j, 8])
                        )  ## strain measured using the HD
                        uncertaintyPts[j, 8] = ustrain(
                            (
                                (pCstInkeVS * speedLightInAPerS)
                                / (
                                    2
                                    * d0[peakNumber]
                                    * np.sin(np.deg2rad(0.5 * AngleHD))
                                )
                            ),
                            np.polyval(calibCoeffsHD, allPtsInPeak[j, 8]),
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
                        and np.polyval(calibCoeffsVD, allPtsInPeak[j, 8]) > 0
                    ):  ## it is the vertical detector
                        pts[j, 8] = np.log(
                            (
                                (pCstInkeVS * speedLightInAPerS)
                                / (
                                    2
                                    * d0[peakNumber]
                                    * np.sin(np.deg2rad(0.5 * AngleVD))
                                )
                            )
                            / np.polyval(calibCoeffsVD, allPtsInPeak[j, 8])
                        )  ## strain measured using the VD
                        uncertaintyPts[j, 8] = ustrain(
                            (
                                (pCstInkeVS * speedLightInAPerS)
                                / (
                                    2
                                    * d0[peakNumber]
                                    * np.sin(np.deg2rad(0.5 * AngleVD))
                                )
                            ),
                            np.polyval(calibCoeffsVD, allPtsInPeak[j, 8]),
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
                    qq1, qq2, qq3 = diffVector(allPtsInPeak[j, 3:8])
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

    h5Save.close()


if __name__ == "__main__":
    run_from_cli(preStraind0cstEDD)
