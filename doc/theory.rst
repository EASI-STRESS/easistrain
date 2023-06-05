Theory
======

Scattering vector
-----------------

The wave vector :math:`\vec{k}` of a beam indicates the beam direction. Its magnitude is defined
to be proportional to the spatial frequency (inverse of the wavelength :math:`\lambda`)

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

The scattered intensity from a crystal is maximal in the directions where the lattice planes fulfill the 
Laue conditions for diffraction

.. math::

    \vec{q} = 2\pi m \vec{n}_H \quad m \in \mathbb{Z}

When only considering the vector lengths you get Bragg's law

.. math::

    2d\sin\theta = m\lambda = m\frac{hc}{E}

where :math:`\lambda` and :math:`E` the incident beam wavelength and energy respectively.

Strain tensor
-------------

The strain tensor in the direction of the scattering vector is

.. math::

    \epsilon^{hkl}(\phi, \psi) = \ln\left(\frac{d^{hkl}(\phi, \psi)}{d^{hkl}_0(\phi, \psi)}\right)

where :math:`\phi` and :math:`\psi` are the azimuth and polar angle of the scattering vector :math:`\vec{q}^{hkl}` in the sample reference frame.

For angular dispersive diffraction (fix energy)

.. math::

    \epsilon^{hkl}(\phi, \psi) = \ln\left(\frac{\sin\theta^{hkl}_0(\phi, \psi)}{\sin\theta^{hkl}(\phi, \psi)}\right)

For energy dispersive diffraction (fix angle)

.. math::

    \epsilon^{hkl}(\phi, \psi) = \ln\left(\frac{E^{hkl}_0(\phi, \psi)}{E^{hkl}(\phi, \psi)}\right)
