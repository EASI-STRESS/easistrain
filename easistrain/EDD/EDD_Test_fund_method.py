import numpy as np
import scipy.optimize



########### Example of the arguments of the main function #############
def diffVector(angles, e11, e22, e33, e23, e13, e12):
    phi = np.radians(angles[:, 0])
    chi = np.radians(angles[:, 1])
    omega = np.radians(angles[:, 2])
    theta = np.radians(angles[:, 3])
    delta = np.radians(angles[:, 4])
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
    defDirMeas = (
        (e11 * q1 ** 2)
        + (e22 * q2 ** 2)
        + (e33 * q3 ** 2)
        + (2 * e12 * q1 * q2)
        + (2 * e13 * q1 * q3)
        + (2 * e23 * q2 * q3)
    )
    return q1, q2, q3, defDirMeas


def deforDirMeas(angles, e11, e22, e33, e23, e13, e12):
    phi = np.radians(angles[:, 0])
    chi = np.radians(angles[:, 1])
    omega = np.radians(angles[:, 2])
    theta = np.radians(angles[:, 3])
    delta = np.radians(angles[:, 4])
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
    defDirMeas = (
        (e11 * q1 ** 2)
        + (e22 * q2 ** 2)
        + (e33 * q3 ** 2)
        + (2 * e12 * q1 * q2)
        + (2 * e13 * q1 * q3)
        + (2 * e23 * q2 * q3)
    )
    return defDirMeas
