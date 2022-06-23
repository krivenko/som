.. _observables:

Supported observables
=====================

.. highlight:: none
.. currentmodule:: som

This page lists all kinds of dynamical observables currently supported by SOM
and explicitly states the integral equation being solved for each of them.

Thermal Green's function of fermions
------------------------------------

.. _fermiongf:

A Green's function of fermions at temperature :math:`T = 1/\beta` is defined
as

.. math::

    G(\tau) = -\langle \mathbb{T}_\tau c(\tau) c^\dagger(0)\rangle.

Its real-frequency counterpart is the retarded Green's function
:math:`G^\mathrm{ret}(\epsilon)`, which is directly connected to the spectral
function :math:`A(\epsilon)`,

.. math::

    A(\epsilon) = -\frac{1}{\pi}\Im G^\mathrm{ret}(\epsilon)

and can be recovered from it using the Hilbert transform,

.. math::

    G^\mathrm{ret}(\epsilon) = -\int\limits_{-\infty}^\infty
    d\epsilon' \frac{A(\epsilon')}{\epsilon' - \epsilon - i0}.

The spectral function is normalized according to

.. math::
    \mathcal{N} = \int\limits_{-\infty}^\infty d\epsilon A(\epsilon).

.. note::

    :math:`\mathcal{N} = 1` for a fermionic Green's function.

    Using the same integral equations, one can also continue fermionic
    self-energies as long as they do not contain a static Hartree-Fock
    contribution (i.e. they decay to 0 as :math:`\omega\to\infty`).
    In this case norms must be computed separately by the user as first spectral
    moments of the self-energy. For derivation of the spectral moments see,
    for instance,

    `M. Potthoff, T. Wegner, and W. Nolting, Phys. Rev. B 55, 16132 (1997) <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.55.16132>`_.

This observable kind is selected and one of the following integral equations
is solved when the :class:`Som` object is constructed with ``kind="FermionGf"``.

.. list-table::
    :header-rows: 1
    :widths: 25 25 50

    * - Imaginary time
      - Imaginary frequencies
      - Legendre orthogonal polynomials
    * -
        .. math::
            G(\tau) = -\int\limits_{-\infty}^\infty d\epsilon
                \frac{e^{-\tau\epsilon}}{1+e^{-\beta\epsilon}}
                A(\epsilon).
      -
        .. math::
            G(i\omega_n) = \int\limits_{-\infty}^\infty d\epsilon
                \frac{1}{i\omega_n-\epsilon}
                A(\epsilon).
      -
        .. math::
            G(\ell) = -\int\limits_{-\infty}^\infty d\epsilon
                \frac{\beta\sqrt{2\ell+1} (-\mathrm{sgn}(\epsilon))^\ell
                    i_{\ell}(\beta|\epsilon|/2)}
                {2\cosh(\beta\epsilon/2)}
                A(\epsilon),

        where :math:`i_\ell(x)` is the modified spherical Bessel function of
        the first kind.

.. _fermiongfsymm:

In many model calculations fermionic Green's functions obey the particle-hole
symmetry, which manifests itself in the symmetry of the spectral  function,
:math:`A(-\epsilon) = A(\epsilon)`.
One can enforce this symmetry by constructing :class:`Som` with
``kind="FermionGfSymm"``. More precisely, SOM will

- Use symmetrized integral kernels (see below) sensitive only to the symmetric
  part of :math:`A(\epsilon)`;
- Symmetrize calculated :ref:`final solution <final_solution>`
  :math:`A(\epsilon)` before recovering :math:`G^\mathrm{ret}(\epsilon)`.

.. list-table::
    :header-rows: 1
    :widths: 25 25 50

    * - Imaginary time
      - Imaginary frequencies
      - Legendre orthogonal polynomials
    * -
        .. math::
            G(\tau) = -\int\limits_{-\infty}^\infty \frac{d\epsilon}{2}
                \frac{e^{-\tau\epsilon} + e^{-(\beta-\tau)\epsilon}}
                {1+e^{-\beta\epsilon}}
                A(\epsilon).
      -
        .. math::
            G(i\omega_n) = -\int\limits_{-\infty}^\infty d\epsilon
                \frac{i\omega_n}{\omega_n^2+\epsilon^2}
                A(\epsilon).
      -
        .. math::
            G(\ell) = \left\{
                \begin{array}{ll}
                -\int\limits_{-\infty}^\infty
                d\epsilon
                \frac{\beta\sqrt{2\ell+1} i_{\ell}(\beta|\epsilon|/2)}
                {2\cosh(\beta\epsilon/2)} A(\epsilon), &\ell\ \mathrm{ even},\\
                0, &\ell\ \mathrm{odd}.
                \end{array}\right.,

        where :math:`i_\ell(x)` is the modified spherical Bessel function of
        the first kind.

Thermal Green's function of bosons, dynamical susceptibilities and conductivity
-------------------------------------------------------------------------------

.. _bosoncorr:

The following correlators share the same general form
:math:`\chi_{OO^\dagger}(\tau) =
\langle \mathbb{T}_\tau \hat O(\tau) \hat O^\dagger(0)\rangle`,
where :math:`\hat O` is a bosonic or boson-like (fermion-number-conserving)
operator.

- Matsubara Green's function of bosons :math:`G(\tau) =
  \langle \mathbb{T}_\tau a(\tau) a^\dagger(0)\rangle`;
- Charge susceptibility :math:`\chi_{NN}(\tau) =
  \langle \mathbb{T}_\tau \hat N(\tau) \hat N(0)\rangle`;
- Longitudinal magnetic susceptibility :math:`\chi_{zz}(\tau) =
  \langle \mathbb{T}_\tau \hat S_z(\tau) \hat S_z(0)\rangle`;
- Transverse magnetic susceptibility :math:`\chi_{-+}(\tau) =
  \langle \mathbb{T}_\tau \hat S_-(\tau) \hat S_+(0)\rangle`;
- Optical conductivity :math:`\sigma(\tau) =
  \langle \mathbb{T}_\tau \hat j(\tau) \hat j(0)\rangle`.

The real-time and real-frequency counterparts of :math:`\chi_{OO^\dagger}(\tau)`
are

.. math::

    \chi_{OO^\dagger}(t) =
        i\theta(t)\langle[\hat O(t),\hat O^\dagger(0)]\rangle,

    \chi_{OO^\dagger}(\epsilon) =
        \int\limits_{-\infty}^\infty dt\ e^{i\epsilon t}\chi_{OO^\dagger}(t).

The imaginary part of :math:`\chi_{OO^\dagger}(\epsilon)` obeys
:math:`\mathrm{sgn}(\Im\chi_{OO^\dagger}(\epsilon)) = \mathrm{sgn}(\epsilon)`,
which allows to introduce a non-negative auxiliary function
:math:`A(\epsilon) = \Im\chi_{OO^\dagger}(\epsilon) / \epsilon`. It plays the
role of the spectral function for this class of continuation problem.

Norm of the spectral function is defined as

.. math::
    \mathcal{N} = \int\limits_{-\infty}^\infty d\epsilon A(\epsilon) =
        \pi\chi_{OO^\dagger}(i\Omega_n=0).

The correlator of a real frequency is recovered according to

    .. math::
        \chi_{OO^\dagger}(\epsilon) = \frac{1}{\pi}\int\limits_{-\infty}^\infty
        d\epsilon' \frac{\epsilon' A(\epsilon')}{\epsilon' - \epsilon - i0}.

This observable kind is selected and one of the following integral equations
is solved when the :class:`Som` object is constructed with ``kind="BosonCorr"``.

.. list-table::
    :header-rows: 1
    :widths: 25 25 50

    * - Imaginary time
      - Imaginary frequencies
      - Legendre orthogonal polynomials
    * -
        .. math::
            \chi_{OO^\dagger}(\tau) =
                \int\limits_{-\infty}^\infty \frac{d\epsilon}{\pi}
                \frac{\epsilon e^{-\tau\epsilon}}{1-e^{-\beta\epsilon}}
                A(\epsilon).
      -
        .. math::
            \chi_{OO^\dagger}(i\Omega_n) = \int\limits_{-\infty}^\infty
                \frac{d\epsilon}{\pi}\frac{-\epsilon}{i\Omega_n - \epsilon}
                A(\epsilon).
      -
        .. math::
            \chi_{OO^\dagger}(\ell) = \int\limits_{-\infty}^\infty
                \frac{d\epsilon}{\pi}
                \frac{\beta\epsilon\sqrt{2\ell+1}(-\mathrm{sgn}(\epsilon))^\ell
                i_{\ell}(\beta|\epsilon|/2)}
                {2\sinh(\beta\epsilon/2)}
                A(\epsilon),

        where :math:`i_\ell(x)` is the modified spherical Bessel function of
        the first kind.

.. _bosonautocorr:

When :math:`\hat O^\dagger = \hat O`, such as in the case of the charge
susceptibility, the longitudinal magnetic susceptibility and the optical
conductivity, it is **recommended** to exploit the additional spectral symmetry
:math:`A(-\epsilon) = A(\epsilon)`. If :class:`Som` is constructed with
``kind="BosonAutoCorr"``, it

- Uses symmetrized integral kernels (see below) sensitive only to the symmetric
  part of :math:`A(\epsilon)`;
- Symmetrizes calculated :ref:`final solution <final_solution>`
  :math:`A(\epsilon)` before recovering :math:`\chi_{OO^\dagger}(\epsilon)`.

.. list-table::
    :header-rows: 1
    :widths: 25 25 50

    * - Imaginary time
      - Imaginary frequencies
      - Legendre orthogonal polynomials
    * -
        .. math::
            \chi_{OO}(\tau) = \int\limits_{-\infty}^\infty
                \frac{d\epsilon}{2\pi}
                \frac{\epsilon (e^{-\tau\epsilon}+e^{-(\beta-\tau)\epsilon})}
                {1-e^{-\beta\epsilon}}
                A(\epsilon).

      -
        .. math::
            \chi_{OO}(i\Omega_n) = \int\limits_{-\infty}^\infty
                \frac{d\epsilon}{\pi}
                \frac{\epsilon^2}{\Omega_n^2+\epsilon^2}
                A(\epsilon).
      -
        .. math::
                \chi_{OO}(\ell) = \left\{
                    \begin{array}{ll}
                    \int\limits_{-\infty}^\infty
                    \frac{d\epsilon}{\pi}
                    \frac{\beta\epsilon\sqrt{2\ell+1}
                    i_{\ell}(\beta|\epsilon|/2)}
                    {2\sinh(\beta\epsilon/2)}
                        A(\epsilon),&\ell\ \mathrm{ even},\\
                    0, &\ell\ \mathrm{odd}.
                    \end{array}\right.,

        where :math:`i_\ell(x)` is the modified spherical Bessel function of
        the first kind.

Dynamical response functions at zero temperature
------------------------------------------------

.. _zerotemp:

In the limit of zero temperature (:math:`\beta=1/T\to\infty`), both Fermi-Dirac
and Bose-Einstein distributions approach the shape of a step function
:math:`\theta(\epsilon)`. It is, therefore, sufficient to consider a spectral
function :math:`A(\epsilon)` defined for :math:`\epsilon\geq 0` regardless of
operator statistics. Formally, the Matsubara time segment
:math:`\tau\in[0; \beta]` becomes infinitely long in this limit. Practically,
however, it is still possible to consider :math:`\tau` varying on a reduced
segment :math:`[0; \tau_\mathrm{max}]`. It is also possible
to introduce a countable sequence of fictitious Matsubara frequencies
:math:`\omega_n = (2n+1)\pi/\tau_\mathrm{max}` or
:math:`\omega_n = 2n\pi/\tau_\mathrm{max}`.

In the :math:`T=0` case, SOM defines the non-negative spectral function as
:math:`A(\epsilon) = -(1/\pi)\Im G(\epsilon)` and its norm as

.. math::
    \mathcal{N} = \int\limits_0^\infty d\epsilon A(\epsilon).

The original response function of a real frequency is recovered according to

.. math::
        G(\epsilon) = -\int\limits_0^\infty
        d\epsilon' \frac{A(\epsilon')}{\epsilon' - \epsilon - i0}.

This observable kind is selected and one of the following integral equations
is solved when the :class:`Som` object is constructed with ``kind="ZeroTemp"``.

.. list-table::
    :header-rows: 1
    :widths: 25 25 50

    * - Imaginary time
      - Imaginary frequencies
      - Legendre orthogonal polynomials
    * -
        .. math::
            G(\tau) = -\int\limits_0^\infty d\epsilon
                e^{-\tau\epsilon}
                A(\epsilon).

      -
        .. math::
            G(i\omega_n) = \int\limits_0^\infty d\epsilon
                \frac{1}{i\omega_n-\epsilon}
                A(\epsilon).
      -
        .. math::
            G(\ell) = \int\limits_0^\infty d\epsilon
                \tau_\mathrm{max}(-1)^{\ell+1}\sqrt{2\ell+1}
                i_{\ell}\left(\frac{\epsilon\tau_\mathrm{max}}{2}\right)
                \exp\left(-\frac{\epsilon\tau_\mathrm{max}}{2}\right)
                A(\epsilon),

        where :math:`i_\ell(x)` is the modified spherical Bessel function of
        the first kind.
