import numpy


def diffVector(angles, e11, e22, e33, e23, e13, e12):
    phi = numpy.radians(angles[:, 0])
    chi = numpy.radians(angles[:, 1])
    omega = numpy.radians(angles[:, 2])
    theta = numpy.radians(angles[:, 3])
    delta = numpy.radians(angles[:, 4])
    q1 = (
        (numpy.cos(theta) * numpy.cos(chi) * numpy.sin(delta) * numpy.sin(phi))
        + (
            numpy.cos(delta)
            * numpy.cos(theta)
            * (
                (numpy.cos(phi) * numpy.sin(omega))
                - (numpy.cos(omega) * numpy.sin(phi) * numpy.sin(chi))
            )
        )
        - numpy.sin(theta)
        * (
            (numpy.cos(phi) * numpy.cos(omega))
            + (numpy.sin(phi) * numpy.sin(chi) * numpy.sin(omega))
        )
    )
    q2 = (
        numpy.cos(delta)
        * numpy.cos(theta)
        * (
            (numpy.cos(phi) * numpy.cos(omega) * numpy.sin(chi))
            + (numpy.sin(phi) * numpy.sin(omega))
        )
        - (numpy.cos(theta) * numpy.cos(phi) * numpy.cos(chi) * numpy.sin(delta))
        - (
            numpy.sin(theta)
            * (
                (numpy.cos(omega) * numpy.sin(phi))
                - (numpy.cos(phi) * numpy.sin(chi) * numpy.sin(omega))
            )
        )
    )
    q3 = (
        (numpy.cos(delta) * numpy.cos(theta) * numpy.cos(chi) * numpy.cos(omega))
        + (numpy.cos(theta) * numpy.sin(delta) * numpy.sin(chi))
        + (numpy.cos(chi) * numpy.sin(theta) * numpy.sin(omega))
    )
    defDirMeas = (
        (e11 * q1**2)
        + (e22 * q2**2)
        + (e33 * q3**2)
        + (2 * e12 * q1 * q2)
        + (2 * e13 * q1 * q3)
        + (2 * e23 * q2 * q3)
    )
    return q1, q2, q3, defDirMeas


def deforDirMeas(angles, e11, e22, e33, e23, e13, e12):
    phi = numpy.radians(angles[:, 0])
    chi = numpy.radians(angles[:, 1])
    omega = numpy.radians(angles[:, 2])
    theta = numpy.radians(angles[:, 3])
    delta = numpy.radians(angles[:, 4])
    q1 = (
        (numpy.cos(theta) * numpy.cos(chi) * numpy.sin(delta) * numpy.sin(phi))
        + (
            numpy.cos(delta)
            * numpy.cos(theta)
            * (
                (numpy.cos(phi) * numpy.sin(omega))
                - (numpy.cos(omega) * numpy.sin(phi) * numpy.sin(chi))
            )
        )
        - numpy.sin(theta)
        * (
            (numpy.cos(phi) * numpy.cos(omega))
            + (numpy.sin(phi) * numpy.sin(chi) * numpy.sin(omega))
        )
    )
    q2 = (
        numpy.cos(delta)
        * numpy.cos(theta)
        * (
            (numpy.cos(phi) * numpy.cos(omega) * numpy.sin(chi))
            + (numpy.sin(phi) * numpy.sin(omega))
        )
        - (numpy.cos(theta) * numpy.cos(phi) * numpy.cos(chi) * numpy.sin(delta))
        - (
            numpy.sin(theta)
            * (
                (numpy.cos(omega) * numpy.sin(phi))
                - (numpy.cos(phi) * numpy.sin(chi) * numpy.sin(omega))
            )
        )
    )
    q3 = (
        (numpy.cos(delta) * numpy.cos(theta) * numpy.cos(chi) * numpy.cos(omega))
        + (numpy.cos(theta) * numpy.sin(delta) * numpy.sin(chi))
        + (numpy.cos(chi) * numpy.sin(theta) * numpy.sin(omega))
    )
    defDirMeas = (
        (e11 * q1**2)
        + (e22 * q2**2)
        + (e33 * q3**2)
        + (2 * e12 * q1 * q2)
        + (2 * e13 * q1 * q3)
        + (2 * e23 * q2 * q3)
    )
    return defDirMeas
