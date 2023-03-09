Coordinate systems
==================

This page defines the coordinate systems in which *easistrain* defines the scattering vector and .

Sample reference frame
----------------------

The cartesian coordinates :math:`Q` of the normalized scattering vector :math:`\hat{q}` in the sample reference frame :math:`(S_1, S_2, S_3)` are given by

.. math::

    Q = \begin{bmatrix}
        \cos \phi \sin \psi \\
        \sin \phi \sin \psi \\
        \cos \psi
        \end{bmatrix}

where :math:`\psi` the polar angle and :math:`\phi` the azimuth angle in spherical coordinates.

Scattering vector and diffraction
---------------------------------

The wave vector :math:`\vec{k}` of a beam indicates the beam direction. Its magnitude is defined
to be proportional to the spatial frequence (inverse of the wavelength :math:`\lambda`)

.. math::

    \|\vec{k}\|= \frac{2\pi}{\lambda}

The scattering vector :math:`\vec{q}` is defined as the difference between the wave vector of the
scattered and incident beam

.. math::

   \vec{q} = \vec{k}_s - \vec{k}_0

The scattering angle :math:`2\theta` between incident and scattered beam is

.. math::

    \cos 2\theta = \hat{k}_s\cdot\hat{k}_0

For elastic scattering :math:`\lambda_s = \lambda_0` the length of the scattering vector
is related to the scattering angle as follows

.. math::

    \|\vec{q}\|= \frac{4\pi}{\lambda} \sin\theta

We define the normal vector :math:`\vec{n}_H` to a family of lattice planes :math:`H = \{hkl\}` as

.. math::

   \| \vec{n}_H \|= \frac{1}{d_H}

The scattered intensity from a crytal is maximal in the directions where the lattice planes fullfill the 
Laue conditions for diffraction

.. math::

    \vec{q} = 2\pi m \vec{n}_H \quad m \in \mathbb{Z}

When only considering the vector lengths you get Bragg's law

.. math::

    2d\sin\theta = m\lambda = m\frac{hc}{E}

where :math:`\lambda` and :math:`E` the incident beam wavelength and energy respectively.
