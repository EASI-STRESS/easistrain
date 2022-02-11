from typing import Sequence, Tuple
import numpy as np
import h5py
import scipy.optimize
from easistrain.EDD.math import compute_qs
from easistrain.EDD.utils import run_from_cli


def guess_strain(
    meas_strain, scattering_x, scattering_y, scattering_z
) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """Initial guess of the strain tensor"""
    scattering_yz = scattering_y * scattering_z
    scattering_xz = scattering_x * scattering_z
    scattering_xy = scattering_x * scattering_y

    raw_eps11_guess = meas_strain[scattering_x >= 0.9]
    guessEps11 = np.mean(raw_eps11_guess) if len(raw_eps11_guess) > 0 else 0

    raw_eps22_guess = meas_strain[scattering_y >= 0.9]
    guessEps22 = np.mean(raw_eps22_guess) if len(raw_eps22_guess) > 0 else 0

    raw_eps33_guess = meas_strain[scattering_z >= 0.9]
    guessEps33 = np.mean(raw_eps33_guess) if len(raw_eps33_guess) > 0 else 0

    raw_eps23_guess = meas_strain[scattering_yz >= 0.9]
    guessEps23 = (
        np.mean(raw_eps23_guess) - 0.5 * (guessEps22 + guessEps33)
        if len(raw_eps23_guess) > 0
        else 0
    )

    raw_eps13_guess = meas_strain[scattering_xz >= 0.9]
    guessEps13 = (
        np.mean(raw_eps13_guess) - 0.5 * (guessEps11 + guessEps33)
        if len(raw_eps13_guess) > 0
        else 0
    )  ## guess of strainxz

    raw_eps12_guess = meas_strain[scattering_xy >= 0.9]
    guessEps12 = (
        raw_eps12_guess - 0.5 * (guessEps11 + guessEps22)
        if len(raw_eps12_guess) > 0
        else 0
    )

    max_strain = np.ones((6)) * (10**-10)

    for j, meas_e in enumerate(
        [
            scattering_x,
            scattering_y,
            scattering_z,
            scattering_yz,
            scattering_xz,
            scattering_xy,
        ]
    ):
        if len(meas_strain[meas_e >= 0.1]) > 0:
            max_strain[j] = np.inf

    return (
        np.array(
            [
                guessEps11,
                guessEps22,
                guessEps33,
                guessEps23,
                guessEps13,
                guessEps12,
            ]
        ),
        (-max_strain, max_strain),
    )


def guess_stress(
    strain_guess, strain_bounds, XEC0, XEC1
) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray]]:

    module_guess: np.ndarray = np.zeros_like(strain_guess)
    module_guess[3:] = (-XEC0 / (XEC1 + 3 * XEC0)) * np.sum(strain_guess[3:])

    stress_guess = (1 / XEC1) * (strain_guess + module_guess)

    min_strain, max_strain = strain_bounds
    max_stress: np.ndarray = np.zeros_like(max_strain)
    max_stress[3:] = max_strain[3:]
    max_stress[:3] = np.inf
    return (
        stress_guess,
        (-max_stress, max_stress),
    )


def strain_in_meas_direction(angles, e11, e22, e33, e23, e13, e12):
    q1, q2, q3 = compute_qs(angles)
    return (
        (e11 * q1**2)
        + (e22 * q2**2)
        + (e33 * q3**2)
        + (2 * e12 * q1 * q2)
        + (2 * e13 * q1 * q3)
        + (2 * e23 * q2 * q3)
    )


def stress_in_meas_direction(anglesAndXEC, s11, s22, s33, s23, s13, s12):
    q1, q2, q3 = compute_qs(anglesAndXEC[:-1])
    S1 = anglesAndXEC[-1, 0]
    dS2 = anglesAndXEC[-1, 1]

    return (S1 * (s11 + s22 + s33)) + (
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


def strainStressTensor(
    fileRead: str, fileSave: str, numberOfPeaks: int, XEC: Sequence[float]
):

    if len(XEC) < numberOfPeaks * 2:
        raise ValueError(
            f"XEC must have a length of numberOfPeaks*2 ({numberOfPeaks*2})"
        )

    h5Save = h5py.File(fileSave, "a")  ## create/append h5 file to save in

    with h5py.File(fileRead, "r") as h5Read:
        for peakNumber in range(numberOfPeaks):
            peak_group = h5Save.create_group(
                f"peak_{str(peakNumber).zfill(4)}"
            )  ## create group for each peak in the strainTensor group
            input_peak_group = h5Read[f"STRAIN_with_d0/peak_{str(peakNumber).zfill(4)}"]
            assert isinstance(input_peak_group, h5py.Group)
            for i in range(len(input_peak_group) // 2):
                input_point_data = input_peak_group[f"point_{str(i).zfill(5)}"][()]
                assert isinstance(input_point_data, np.ndarray)
                input_point_errors = input_peak_group[
                    f"uncertainty_point_{str(i).zfill(5)}"
                ][()]
                assert isinstance(input_point_errors, np.ndarray)

                meas_angles = input_point_data[:, 3:8]
                meas_strain = input_point_data[:, 8]
                scattering_x = input_point_data[:, 9]
                scattering_y = input_point_data[:, 10]
                scattering_z = input_point_data[:, 11]

                strain_tensor_guess, strain_tensor_bounds = guess_strain(
                    meas_strain, scattering_x, scattering_y, scattering_z
                )
                uncertainty_meas_strain = input_point_errors[:, 8]
                strain_tensor_fit, covarStrains = scipy.optimize.curve_fit(
                    f=strain_in_meas_direction,
                    xdata=meas_angles,
                    ydata=meas_strain,
                    p0=strain_tensor_guess,
                    sigma=uncertainty_meas_strain,
                    bounds=strain_tensor_bounds,
                )

                stress_tensor_guess, stress_tensor_bounds = guess_stress(
                    strain_tensor_guess,
                    strain_tensor_bounds,
                    XEC[2 * peakNumber],
                    XEC[2 * peakNumber + 1],
                )
                stress_tensor_fit, covarStress = scipy.optimize.curve_fit(
                    f=stress_in_meas_direction,
                    xdata=np.append(
                        meas_angles,
                        [[XEC[2 * peakNumber], XEC[2 * peakNumber + 1], 0, 0, 0]],
                        axis=0,
                    ),
                    ydata=meas_strain,
                    p0=stress_tensor_guess,
                    sigma=uncertainty_meas_strain,
                    bounds=stress_tensor_bounds,
                )

                point_name = f"point_{str(i).zfill(5)}"
                point_in_peak_group = peak_group.create_group(point_name)
                point_in_peak_group.create_dataset(
                    "position", data=input_point_data[0, 0:3]
                )
                point_in_peak_group.create_dataset(
                    "position_errors", data=input_point_errors[0, 0:3]
                )
                point_in_peak_group.create_dataset(
                    "strain_tensor_fit", data=strain_tensor_fit
                )
                point_in_peak_group.create_dataset(
                    "strain_tensor_errors",
                    data=(
                        np.sqrt(np.diag(covarStrains))
                        if np.sum(covarStrains) != np.inf
                        else np.mean(input_point_errors[:, 8])
                    ),
                )

                point_in_peak_group.create_dataset(
                    "stress_tensor_fit", data=stress_tensor_fit
                )
                point_in_peak_group.create_dataset(
                    "stress_tensor_errors",
                    data=(
                        np.sqrt(np.diag(covarStress))
                        if np.sum(covarStress) != np.inf
                        else (1 / (XEC[2 * peakNumber + 1] + XEC[2 * peakNumber]))
                        * np.mean(input_point_errors[:, 8])
                    ),
                )

    h5Save.close()
    return


if __name__ == "__main__":
    run_from_cli(strainStressTensor)
