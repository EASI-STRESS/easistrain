Coordinate systems
==================

This page defines the coordinate systems in which *easistrain* defines the scattering vector.

Sample reference frame
++++++++++++++++++++++

Right-handed Eucledian reference frame :math:`(\hat{S}_1, \hat{S}_2, \hat{S}_3)` with

- :math:`\hat{S}_1`: physically meaningfully direction (rolling direction, machining direction, welding direction, ...)
- :math:`\hat{S}_3`: sample surface normal

Laboratory reference frame
++++++++++++++++++++++++++

Right-handed Eucledian reference frame :math:`(\hat{L}_1, \hat{L}_2, \hat{L}_3)` with

- :math:`\hat{L}_1`: beam direction
- :math:`\hat{L}_3`: inverse direction of gravity

Note that :math:`\hat{L}_2` points to the left when looking downstream and :math:`\hat{L}_3` points upwards.

In the code this reference frame is also referred to as the goniometer reference frame.

Scattering vector
+++++++++++++++++

The cartesian coordinates :math:`Q` of the normalized scattering vector :math:`\hat{q}` in the sample reference frame :math:`(\hat{S}_1, \hat{S}_2, \hat{S}_3)` are given by

.. math::

    Q_\text{sample} = \begin{bmatrix}
        \cos \phi \sin \psi \\
        \sin \phi \sin \psi \\
        \cos \psi
        \end{bmatrix}

- :math:`\phi`: azimuth angle of :math:`\hat{q}` in spherical coordinates
- :math:`\psi`: polar angle of :math:`\hat{q}` in spherical coordinates

The cartesian coordinates :math:`Q` of the normalized scattering vector :math:`\hat{q}` in the laboratory reference frame :math:`(\hat{L}_1, \hat{L}_2, \hat{L}_3)` are given by

.. math::

    Q_\text{lab} = \begin{bmatrix}
        -\sin \theta \\
        -\cos \theta \sin \delta \\
        \cos \theta \cos \delta
        \end{bmatrix}

- :math:`\theta`: half the angle between incident beam (along :math:`\hat{L}_1`) and scattered beam
- :math:`\delta`: angle of scattered beam in the :math:`\hat{L}_3\times-\hat{L}_2` plane (careful: :math:`\delta=0^\circ` is :math:`\hat{L}_3` and
:math:`\delta=90^\circ` is :math:`-\hat{L}_2`)

The *easistrain* project refers to a horizontal and vertical detector as

- horizontal detector (:math:`\delta=-90^\circ`): scattering plane is in the synchrotron plane and the detector is positioned on the left when looking downstream

.. math::

    Q_\text{lab,H} = \begin{bmatrix}
        -\sin \theta_H \\
        \cos \theta_H \\
        0
        \end{bmatrix}

- vertical detector (:math:`\delta=0^\circ`): scattering plane is perpendicular to the synchrotron plane and the detector is positioned above the synchrotron plane

.. math::

    Q_\text{lab,V} = \begin{bmatrix}
        -\sin \theta \\
        0 \\
        \cos \theta
        \end{bmatrix}

The activate transformations of :math:`\hat{q}` between cartesian coordinates in laboratory to sample frame are given by

.. math::

    Q_\text{lab} = R_2(-\omega) \cdot R_1(\chi) \cdot R_3(-\varphi) \cdot Q_\text{sample}

.. math::

    Q_\text{sample} = R_3(\varphi) \cdot R_1(-\chi) \cdot R_2(\omega) \cdot Q_\text{lab}

- :math:`\chi`: activate rotation of the sample around :math:`\hat{L}_1`
- :math:`\omega`: activate rotation of the sample around :math:`-\hat{L}_2`
- :math:`\varphi`: activate rotation of the sample around :math:`-\hat{L}_3`
- :math:`R_i`: activate transformation matrix around axis :math:`i`

For example :math:`\chi=0^\circ`, :math:`\omega=90^\circ` and :math:`\varphi=0^\circ`
positions the sample surface perpendicular to the beam (:math:`\hat{S}_3 = -\hat{L}_1`).
