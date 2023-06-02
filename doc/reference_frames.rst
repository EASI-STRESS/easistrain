Coordinate systems
==================

This page defines the coordinate systems in which *easistrain* defines the scattering vector.

Sample reference frame
++++++++++++++++++++++

Right-handed Eucledian reference frame with

- :math:`\hat{S}_1`: physically meaningfully direction (rolling direction, machining direction, welding direction, ...)
- :math:`\hat{S}_2`:
- :math:`\hat{S}_3`: sample surface normal

Goniometer reference frame
++++++++++++++++++++++++++

Right-handed Eucledian reference frame with

- :math:`\hat{G}_1`: beam direction
- :math:`\hat{G}_2`: parallel to the axial direction of the goniometer (direction of gravity, not sure this is always the case)
- :math:`\hat{G}_3`:

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

The cartesian coordinates :math:`Q` of the normalized scattering vector :math:`\hat{q}` in the goniometer reference frame :math:`(\hat{G}_1, \hat{G}_2, \hat{G}_3)` are given by

.. math::

    Q_\text{gonio} = \begin{bmatrix}
        -\sin \theta \\
        -\cos \theta \sin \delta \\
        \cos \theta \cos \delta
        \end{bmatrix}

- :math:`\theta`: half the angle between incident beam (along :math:`\hat{G}_1`) and scattered beam
- :math:`\delta`: angle of scattered beam in the :math:`\hat{G}_3\times\hat{G}_2` plane looking downstream

The *easistrain* project refers to a horizontal and vertical detector as

- horizontal detector (:math:`\delta=0^\circ`): scattering plane is in the synchrotron plane and the detector is positioned on the left when looking downstream

.. math::

    Q_\text{gonio,H} = \begin{bmatrix}
        -\sin \theta_H \\
        0 \\
        \cos \theta_H
        \end{bmatrix}

- vertical detector (:math:`\delta=-90^\circ`): scattering plane is perpendicular ot the synchrotron plane and the detector is positioned below the synchrotron plane

.. math::

    Q_\text{gonio,V} = \begin{bmatrix}
        -\sin \theta_V \\
        \cos \theta_V \\
        0
        \end{bmatrix}

The activate transformation of :math:`\hat{q}` from goniometer to sample frame is

.. math::

    Q_\text{sample} = R_2(-\omega)\cdot R_1(\chi)\cdot R_3(-\varphi)\cdot Q_\text{gonio}

- :math:`\chi`: rotation around :math:`\hat{G}_1`
- :math:`\omega`: rotation around :math:`-\hat{G}_2`
- :math:`\varphi`: rotation around :math:`-\hat{G}_3`
- :math:`R_i`: activate transformation matrix around axis :math:`\hat{G}_i`

For example :math:`\chi=0^\circ`, :math:`\omega=90^\circ` and :math:`\chi=0^\circ`
positions the sample surface perpendicular to the beam (:math:`\hat{G}_3 = -\hat{S}_1`).
