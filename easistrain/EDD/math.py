import numpy


def angles_to_rad(angles: numpy.ndarray):
    rad_angles = numpy.radians(angles)
    # Convert 2*theta in theta
    rad_angles[..., 4] = 0.5 * rad_angles[..., 4]

    # phi, chi, omega, delta, 2*theta as columns
    return numpy.transpose(rad_angles)


def compute_qs(angles: numpy.ndarray):
    rad_angles = angles_to_rad(angles)
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
