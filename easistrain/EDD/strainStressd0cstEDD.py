from typing import Sequence
import numpy as np
import h5py
import scipy.optimize
from easistrain.EDD.utils import run_from_cli


def angles_to_rad(angles: np.ndarray):
    rad_angles = np.radians(angles)
    # Convert 2*theta in theta
    rad_angles[:, 4] = 0.5 * rad_angles[:, 4]

    # phi, chi, omega, delta, 2*theta as columns
    return np.transpose(rad_angles)


def compute_q_factors(angles: np.ndarray):
    rad_angles = angles_to_rad(angles)
    cos_phi, cos_chi, cos_omega, cos_delta, cos_theta = np.cos(rad_angles)
    sin_phi, sin_chi, sin_omega, sin_delta, sin_theta = np.sin(rad_angles)
    q1 = (
        (cos_theta * cos_chi * sin_delta * sin_phi)
        + (
            cos_delta
            * cos_theta
            * ((cos_phi * sin_omega) - (cos_omega * sin_phi * sin_chi))
        )
        - sin_theta * ((cos_phi * cos_omega) + (sin_phi * sin_chi * sin_omega))
    )
    q2 = (
        cos_delta
        * cos_theta
        * ((cos_phi * cos_omega * sin_chi) + (sin_phi * sin_omega))
        - (cos_theta * cos_phi * cos_chi * sin_delta)
        - (sin_theta * ((cos_omega * sin_phi) - (cos_phi * sin_chi * sin_omega)))
    )
    q3 = (
        (cos_delta * cos_theta * cos_chi * cos_omega)
        + (cos_theta * sin_delta * sin_chi)
        + (cos_chi * sin_theta * sin_omega)
    )
    return q1, q2, q3


def deforDirMeas(angles, e11, e22, e33, e23, e13, e12):
    q1, q2, q3 = compute_q_factors(angles)
    defDirMeas = (
        (e11 * q1**2)
        + (e22 * q2**2)
        + (e33 * q3**2)
        + (2 * e12 * q1 * q2)
        + (2 * e13 * q1 * q3)
        + (2 * e23 * q2 * q3)
    )
    return defDirMeas


def deforDirMeasStress(anglesAndXEC, s11, s22, s33, s23, s13, s12):
    q1, q2, q3 = compute_q_factors(anglesAndXEC[:-1])
    S1 = anglesAndXEC[-1, 0]
    dS2 = anglesAndXEC[-1, 1]
    defDirMeasStress = (S1 * (s11 + s22 + s33)) + (
        dS2
        * (
            (s11 * q1**2)
            + (s22 * q2**2)
            + (s33 * q3**2)
            + (2 * s12 * q1 * q2)
            + (2 * s13 * q1 * q3)
            + (2 * s23 * q2 * q3)
        )
    )
    return defDirMeasStress


def strainStressTensor(
    fileRead: str, fileSave: str, numberOfPeaks: int, XEC: Sequence[float]
):

    if len(XEC) < numberOfPeaks * 2:
        raise ValueError(
            f"XEC must have a length of numberOfPeaks*2 ({numberOfPeaks*2})"
        )

    h5Save = h5py.File(fileSave, "a")  ## create/append h5 file to save in
    strainTensor = h5Save.create_group(
        "strain_tensor"
    )  ## creation of the strain tensor group
    stressTensor = h5Save.create_group(
        "stress_tensor"
    )  ## creation of the stress tensor group

    with h5py.File(fileRead, "r") as h5Read:
        for peakNumber in range(numberOfPeaks):
            peakGroup = strainTensor.create_group(
                f"peak_{str(peakNumber).zfill(4)}"
            )  ## create group for each peak in the strainTensor group
            peakGroupStress = stressTensor.create_group(
                f"peak_{str(peakNumber).zfill(4)}"
            )  ## create group for each peak in the stressTensor group
            input_peak_group = h5Read[f"STRAIN_with_d0/peak_{str(peakNumber).zfill(4)}"]
            assert isinstance(input_peak_group, h5py.Group)
            for i in range(len(input_peak_group) // 2):
                pointInPeak = input_peak_group[f"point_{str(i).zfill(5)}"][()]
                assert isinstance(pointInPeak, np.ndarray)
                upointInPeak = input_peak_group[f"uncertainty_point_{str(i).zfill(5)}"][
                    ()
                ]
                assert isinstance(upointInPeak, np.ndarray)

                meas_strain = pointInPeak[:, 8]
                meas_exx = pointInPeak[:, 9]
                meas_eyy = pointInPeak[:, 10]
                meas_ezz = pointInPeak[:, 11]
                meas_eyz = meas_eyy * meas_ezz
                meas_exz = meas_exx * meas_ezz
                meas_exy = meas_exx * meas_eyy
                meas_strain_tensor = np.array(
                    [meas_exx, meas_eyy, meas_ezz, meas_eyz, meas_exz, meas_exy]
                )

                raw_eps11_guess = meas_strain[meas_exx >= 0.9]
                guessEps11 = np.mean(raw_eps11_guess) if len(raw_eps11_guess) > 0 else 0

                raw_eps22_guess = meas_strain[meas_eyy >= 0.9]
                guessEps22 = np.mean(raw_eps22_guess) if len(raw_eps22_guess) > 0 else 0

                raw_eps33_guess = meas_strain[meas_ezz >= 0.9]
                guessEps33 = np.mean(raw_eps33_guess) if len(raw_eps33_guess) > 0 else 0

                raw_eps23_guess = meas_strain[meas_eyz >= 0.9]
                guessEps23 = (
                    np.mean(raw_eps23_guess) - 0.5 * (guessEps22 + guessEps33)
                    if len(raw_eps23_guess) > 0
                    else 0
                )

                raw_eps13_guess = meas_strain[meas_exz >= 0.9]
                guessEps13 = (
                    np.mean(raw_eps13_guess) - 0.5 * (guessEps11 + guessEps33)
                    if len(raw_eps13_guess) > 0
                    else 0
                )  ## guess of strainxz

                raw_eps12_guess = meas_strain[meas_exy >= 0.9]
                guessEps12 = (
                    raw_eps12_guess - 0.5 * (guessEps11 + guessEps22)
                    if len(raw_eps12_guess) > 0
                    else 0
                )

                guessSig11 = (1 / XEC[2 * peakNumber + 1]) * (
                    guessEps11
                    + (
                        (
                            -XEC[2 * peakNumber]
                            / (XEC[2 * peakNumber + 1] + 3 * XEC[2 * peakNumber])
                        )
                        * (guessEps11 + guessEps22 + guessEps33)
                    )
                )  ## guess of sigmaxx
                guessSig22 = (1 / XEC[2 * peakNumber + 1]) * (
                    guessEps22
                    + (
                        (
                            -XEC[2 * peakNumber]
                            / (XEC[2 * peakNumber + 1] + 3 * XEC[2 * peakNumber])
                        )
                        * (guessEps11 + guessEps22 + guessEps33)
                    )
                )  ## guess of sigmayy
                guessSig33 = (1 / XEC[2 * peakNumber + 1]) * (
                    guessEps33
                    + (
                        (
                            -XEC[2 * peakNumber]
                            / (XEC[2 * peakNumber + 1] + 3 * XEC[2 * peakNumber])
                        )
                        * (guessEps11 + guessEps22 + guessEps33)
                    )
                )  ## guess of sigmazz
                guessSig23 = (
                    1 / XEC[2 * peakNumber + 1]
                ) * guessEps23  ## guess of sigmayz
                guessSig13 = (
                    1 / XEC[2 * peakNumber + 1]
                ) * guessEps13  ## guess of sigmaxz
                guessSig12 = (
                    1 / XEC[2 * peakNumber + 1]
                ) * guessEps12  ## guess of sigmaxy
                strainTensorInitialGuess = [
                    guessEps11,
                    guessEps22,
                    guessEps33,
                    guessEps23,
                    guessEps13,
                    guessEps12,
                ]  ## initial guess of the strain tensor
                stressTensorInitialGuess = [
                    guessSig11,
                    guessSig22,
                    guessSig33,
                    guessSig23,
                    guessSig13,
                    guessSig12,
                ]  ## initial guess of the stress tensor

                bounds_eps = np.ones((6)) * (10**-10)
                bounds_sig = np.ones((6)) * (10**-10)

                for j, meas_e in enumerate(meas_strain_tensor):
                    if len(meas_strain[meas_e >= 0.1]) > 0:
                        bounds_eps[j] = np.inf
                        bounds_sig[j] = np.inf

                # bounds on sigma11, sigma22, sigma33
                bounds_sig[:3] = np.inf

                strainTensorComponents, covarStrains = scipy.optimize.curve_fit(
                    f=deforDirMeas,
                    xdata=pointInPeak[:, 3:8],
                    ydata=pointInPeak[:, 8],
                    p0=strainTensorInitialGuess,
                    sigma=upointInPeak[:, 8],
                    bounds=(-bounds_eps, bounds_eps),
                )  ## fit of the strain tensor
                stressTensorComponents, covarStress = scipy.optimize.curve_fit(
                    f=deforDirMeasStress,
                    xdata=np.append(
                        pointInPeak[:, 3:8],
                        np.array(
                            [[XEC[2 * peakNumber], XEC[2 * peakNumber + 1], 0, 0, 0]]
                        ),
                        axis=0,
                    ),
                    ydata=pointInPeak[:, 8],
                    p0=stressTensorInitialGuess,
                    sigma=upointInPeak[:, 8],
                    bounds=(-bounds_sig, bounds_sig),
                )  ## fit of the stress tensor

                pointInPeakGroup = peakGroup.create_dataset(
                    f"point_{str(i).zfill(5)}",
                    dtype="float64",
                    data=np.zeros((1, 9), "float64"),
                )  ## strain tensor in point i
                pointInPeakGroupStress = peakGroupStress.create_dataset(
                    f"point_{str(i).zfill(5)}",
                    dtype="float64",
                    data=np.zeros((1, 9), "float64"),
                )  ## stress tensor in point i
                upointInPeakGroup = peakGroup.create_dataset(
                    f"uncertainty_point_{str(i).zfill(5)}",
                    dtype="float64",
                    data=np.zeros((1, 9), "float64"),
                )  ## uncertinty on the strain tensor component in point i
                upointInPeakGroupStress = peakGroupStress.create_dataset(
                    f"uncertainty_point_{str(i).zfill(5)}",
                    dtype="float64",
                    data=np.zeros((1, 9), "float64"),
                )  ## uncertainty on the stress tensor component in point i

                pointInPeakGroup[0, 0:3] = pointInPeak[0, 0:3]
                pointInPeakGroup[0, 3:9] = strainTensorComponents

                upointInPeakGroup[0, 0:3] = upointInPeak[0, 0:3]
                upointInPeakGroup[0, 3:9] = (
                    np.sqrt(np.diag(covarStrains))
                    if np.sum(covarStrains) != np.inf
                    else np.mean(upointInPeak[:, 8])
                )

                pointInPeakGroupStress[0, 0:3] = pointInPeak[0, 0:3]
                pointInPeakGroupStress[0, 3:9] = stressTensorComponents

                upointInPeakGroupStress[0, 0:3] = upointInPeak[0, 0:3]
                upointInPeakGroupStress[0, 3:9] = (
                    np.sqrt(np.diag(covarStress))
                    if np.sum(covarStress) != np.inf
                    else (1 / (XEC[2 * peakNumber + 1] + XEC[2 * peakNumber]))
                    * np.mean(upointInPeak[:, 8])
                )

    h5Save.close()
    return


if __name__ == "__main__":
    run_from_cli(strainStressTensor)
