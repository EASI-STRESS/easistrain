from typing import List
import numpy


def compute_qs(angles: numpy.ndarray):
    rad_angles = numpy.radians(angles).T
    rad_angles[4] /= 2  # convert 2*theta to theta

    cos_phi, cos_chi, cos_omega, cos_delta, cos_theta = numpy.cos(rad_angles)
    sin_phi, sin_chi, sin_omega, sin_delta, sin_theta = numpy.sin(rad_angles)
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


def qsample(angles: numpy.ndarray) -> numpy.ndarray:
    """Normalized scattering vector in the sample frame from
    spherical coordinates in goniometer frame and goniometer rotations.

    Non-optimized implementation of `compute_qs`.
    """
    phi, chi, omega, delta, ttheta = angles
    # Transformation from goniometer to sample frame
    M = roty(-omega).dot(rotx(chi).dot(rotz(-phi)))
    # Coordinates in goniometer frame (spherical to cartesian)
    Q = qgonio(delta, ttheta)
    # Coordinates in ample frame
    return M.T.dot(Q)


def qgonio(delta: float, ttheta: float) -> List[float]:
    """Normalized scattering vector in the goniometer frame: spherical coordinates to cartesian coordinates"""
    rad_delta = numpy.radians(delta)
    rad_theta = numpy.radians(ttheta / 2)
    return [
        -numpy.sin(rad_theta),
        -numpy.cos(rad_theta) * numpy.sin(rad_delta),
        numpy.cos(rad_theta) * numpy.cos(rad_delta),
    ]


def rotx(rx: float) -> numpy.ndarray:
    a = numpy.radians(rx)
    cosa = numpy.cos(a)
    sina = numpy.sin(a)
    return numpy.array(
        [
            [1, 0, 0],
            [0, cosa, -sina],
            [0, sina, cosa],
        ]
    )


def roty(ry: float) -> numpy.ndarray:
    a = numpy.radians(ry)
    cosa = numpy.cos(a)
    sina = numpy.sin(a)
    return numpy.array(
        [
            [cosa, 0, sina],
            [0, 1, 0],
            [-sina, 0, cosa],
        ]
    )


def rotz(rz: float) -> numpy.ndarray:
    a = numpy.radians(rz)
    cosa = numpy.cos(a)
    sina = numpy.sin(a)
    return numpy.array(
        [
            [cosa, -sina, 0],
            [sina, cosa, 0],
            [0, 0, 1],
        ]
    )
