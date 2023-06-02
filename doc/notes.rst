Unordered notes
===============

Task 1: `calibrationEDD`
++++++++++++++++++++++++

Purpose: calibration of the conversion from channels to energy with X-ray source

.. math::

    \text{E} = p_0 x^2 + p_1 x + p_2

Procedure: fit peaks in channels and use their known energy (e.g. Ba or Am fluorescence source) to calculate coefficients :math:`p_i`.

Data:

    - /detectorCalibration/fit_0001_2_1/calibCoeffs/calibCoeffsHD: :math:`p_i` stored in `numpy.polyval` order
    - /detectorCalibration/fit_0001_2_1/calibCoeffs/calibCoeffsVD: :math:`p_i` stored in `numpy.polyval` order
    - /detectorCalibration/fit_0001_2_1/fitLine_{i}/horizontal: NXdata plot of fit
    - /detectorCalibration/fit_0001_2_1/fitLine_{i}/vertical: NXdata plot of fit
    - /detectorCalibration/fit_0001_2_1/fitParams/fitParamsHD: fit parameters
            - columns are ???
            - rows are ???
    - /detectorCalibration/fit_0001_2_1/fitParams/fitParamsVD: fit parameters
            - columns are ???
            - rows are ???

Task 2: `angleCalibEDD`
+++++++++++++++++++++++

Purpose: calibration of the diffraction angle with known sample

.. math::

    2d\sin\theta = \frac{hc}{E}

Procedure: fit peaks in channels and use the energy calibration and known d-spacings to calculate :math:`2\theta`.

Data:

    - /angleCalibration/fit_0001_3/calibratedAngle/calibratedAngleHD: :math:`2\theta` of the horizontal detector
    - /angleCalibration/fit_0001_3/calibratedAngle/calibratedAngleVD: :math:`2\theta` of the vertical detector
    - /angleCalibration/fit_0001_3/fitLine_{}:  NXdata plot of fit
    - /angleCalibration/fit_0001_3/fitParams/fitParamsHD: fit parameters
            - columns are ???
            - rows are ???
    - /angleCalibration/fit_0001_3/fitParams/fitParamsVD: fit parameters
            - columns are ???
            - rows are ???

Task 3: `fitEDD`
++++++++++++++++

Purpose: peak positions in channels for unknown sample

Procedure: fit peaks in channels

Data:

    - /BAIII_AB_4_1_0001_42.1/fit/0000/fitLine_0000: NXdata plot of fit
    - /BAIII_AB_4_1_0001_42.1/fit/0000/fitParams/fitParamsHD: fit parameters
            - columns are ???
            - rows are ???
    - /BAIII_AB_4_1_0001_42.1/fit/0000/fitParams/fitParamsVD: fit parameters
            - columns are ???
            - rows are ???
    - /BAIII_AB_4_1_0001_42.1/positioners: ex, ey, ez, ephi, echi, sy
    - /BAIII_AB_4_1_0001_42.1/tthPositionsGroup/peak_{i}: combine positioners and fit prameters
            - columns are [ex, ey, ez, ephi, echi, sy, delta (azimuth), theta, position (in channel), Intenstity, FWHM, shape factor, goodness factor]
            - 2 rows (horizontal and vertical)

delta is -90 degrees for the horizontal detector and 0 degrees for the vertical detector.

Task 4: `coordTransformation`
+++++++++++++++++++++++++++++

transformation of the coordinates from the gonio reference to the sample reference

Task 5: `regroupPoints`
+++++++++++++++++++++++

regroup all the points


Task 6: `preStraind0cstEDD`
+++++++++++++++++++++++++++

calculates the strain in the measurement direction

Task 7: `strainStressd0cstEDD`
++++++++++++++++++++++++++++++

calculates the strain and stress tensors

