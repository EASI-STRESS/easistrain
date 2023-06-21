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
    """Normalized scattering vector in the sample frame.
    As input we expect the spherical coordinates of the normalized scattering vector
    in the laboratory frame and the orientation of the sample in the laboratory frame.

    Non-optimized implementation of `compute_qs`.
    """
    phi, chi, omega, delta, ttheta = angles
    # Transformation from laboratory to sample frame
    # M = (roty(-omega).dot(rotx(chi)).dot(rotz(-phi))).T
    M = rotz(phi).dot(rotx(-chi)).dot(roty(omega))
    # Coordinates in laboratory frame (spherical to cartesian)
    Q = qlab(delta, ttheta)
    # Coordinates in sample frame
    return M.dot(Q)


def qlab(delta: float, ttheta: float) -> List[float]:
    """Normalized scattering vector in the laboratory frame: spherical coordinates to cartesian coordinates"""
    rad_delta = numpy.radians(delta)
    rad_theta = numpy.radians(ttheta / 2)
    return [
        -numpy.sin(rad_theta),
        -numpy.cos(rad_theta) * numpy.sin(rad_delta),
        numpy.cos(rad_theta) * numpy.cos(rad_delta),
    ]


def rotx(rx: float) -> numpy.ndarray:
    """Active transformation: rotation around X"""
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
    """Active transformation: rotation around Y"""
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
    """Active transformation: rotation around Z"""
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
