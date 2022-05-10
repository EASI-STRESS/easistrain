import scipy.constants

pCstInkeVS = (
    0.001 * scipy.constants.physical_constants["Planck constant in eV s"][0]
)  ## Planck Constant in keV s
speedLightInAPerS = (10 ** 10) * scipy.constants.physical_constants[
    "speed of light in vacuum"
][
    0
]  ## speed of light in Ang/second
