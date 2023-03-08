import numpy


def func_PearsonVII(P, x):
    """PearsonVII function used to describe an X-ray diffraction peak

    - P[0]: Intensity of the peak (maximum intensity of the peak)
    - P[1]: 2Theta (2 * bragg angle)
    - P[2]: FWHM
    - P[3]: Shape factor (Lorentzian: shape factor --> 1 , Gaussian: shape factor --> infinity)
    """
    return (
        P[0]
        / (1 + (2 * (x - P[1]) * (2 ** (1 / P[3]) - 1) ** 0.5) ** 2 / (P[2] ** 2))
        ** P[3]
    )


def func_Gauss(P, x):
    """Gaussian function used to describe an X-ray diffraction peak

    - P[0]: Intensity of the peak (maximum intensity of the peak)
    - P[1]: 2Theta (2 * bragg angle)
    - P[2]: Integral breadth which is related to the FWHM: IB = 0.5*FWHM*sqrt(pi/Ln2) = 1.0644*FWHM
    """
    return P[0] * numpy.exp(-(numpy.pi) * ((x - P[1]) / P[2]) ** 2)


def func_Lorentz(P, x):
    """Lorentz function used to describe an X-ray diffraction peak

    - P[0]: Intensity of the peak (maximum intensity of the peak)
    - P[1]: 2Theta (2 * bragg angle)
    - P[2]: Integral breadth which is related to the FWHM: IB = pi*FWHM/2 = 1.571*FWHM
    """
    return P[0] / (1 + (((numpy.pi) ** 2) * (((x - P[1]) / P[2]) ** 2)))


def func_pseudo_voigt(P, x):
    """Pseudo-Voigt function used to describe an X-ray diffraction peak

    - P[0]: Intensity of the peak (maximum intensity of the peak)
    - P[1]: 2Theta (2 * bragg angle)
    - P[2]: Integral breadth which is related to the FWHM: IB = (0.5*pi*FWHM)/(P[3]+((1-P[3])*sqrt(pi*Ln2)))
    - P[3]: Shape factor (Lorentzian: shape factor --> 1 , Gaussian: shape factor --> 0)
    """
    return P[3] * func_Lorentz(P, x) + ((1 - P[3]) * func_Gauss(P, x))
