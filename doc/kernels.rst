.. _kernels:

Supported integral kernels
==========================

.. note::

    :math:`\beta` denotes the inverse temperature in all equations below.

.. _modified spherical Bessel function of the first kind:
    https://mathworld.wolfram.com/ModifiedSphericalBesselFunctionoftheFirstKind.html

.. _fermiongf:

Fermionic thermal Green's function
----------------------------------

:math:`A(\epsilon) = -(1/\pi)\Im G^\mathrm{ret}(\epsilon)` is the spectral function to be found.

Enabled when `Som()` object is constructed with `kind = "FermionGf"`.

- In imaginary time, :math:`G(\tau)`

    .. math::
        G(\tau) = -\int\limits_{-\infty}^\infty
        d\epsilon \frac{e^{-\tau\epsilon}}{1+e^{-\beta\epsilon}} A(\epsilon).

- At Matsubara frequencies, :math:`G(i\omega_n)`

    .. math::
        G(i\omega_n) = \int\limits_{-\infty}^\infty
        d\epsilon \frac{1}{i\omega_n-\epsilon} A(\epsilon).

- In Legendre polynomial basis, :math:`G(\ell)`

    .. math::
        G(\ell) = -\int\limits_{-\infty}^\infty
        d\epsilon \frac{\beta\sqrt{2\ell+1}(-\mathrm{sgn}(\epsilon))^\ell i_{\ell}(\beta|\epsilon|/2)}
        {2\cosh(\beta\epsilon/2)} A(\epsilon),

  where :math:`i_\ell(x)` is the `modified spherical Bessel function of the first kind`_.

Norm of a solution :math:`A(\epsilon)` is defined as

    .. math::
        N = \int\limits_{-\infty}^\infty d\epsilon A(\epsilon).

The default value :math:`N=1` can be overridden to perform analytical continuation of related fermionic
quantities, such as self-energy.

The retarded Green's function of a real frequency is reconstructed according to

    .. math::
        G^\mathrm{ret}(\epsilon) = -\int\limits_{-\infty}^\infty
        d\epsilon' \frac{A(\epsilon')}{\epsilon' - \epsilon - i0}.

.. _bosoncorr:

Correlator of boson-like operators
----------------------------------

Enabled when `Som()` object is constructed with `kind = "BosonCorr"`.

:math:`A(\epsilon)` is the spectral function to be found. It is defined differently from the fermionic
case, namely :math:`A(\epsilon) = \Im\chi(\epsilon)/\epsilon`, where

    .. math::
        \chi(\epsilon) = \int\limits_{-\infty}^\infty dt\ e^{i\epsilon t}\chi(t),\quad
        \chi(t) = i\theta(t)\langle[\hat O(t),\hat O^\dagger(0)]\rangle.

- In imaginary time, :math:`\chi(\tau)`

    .. math::
        \chi(\tau) = \int\limits_{-\infty}^\infty \frac{d\epsilon}{\pi}
        \frac{\epsilon e^{-\tau\epsilon}}{1-e^{-\beta\epsilon}} A(\epsilon).

- At Matsubara frequencies, :math:`\chi(i\Omega_n)`

    .. math::
        \chi(i\Omega_n) = \int\limits_{-\infty}^\infty
        d\epsilon \frac{1}{\pi}\frac{-\epsilon}{i\Omega_n - \epsilon} A(\epsilon).

- In Legendre polynomial basis, :math:`\chi(\ell)`

    .. math::
        \chi(\ell) = \int\limits_{-\infty}^\infty
        \frac{d\epsilon}{\pi}
        \frac{\beta\epsilon\sqrt{2\ell+1}(-\mathrm{sgn}(\epsilon))^\ell i_{\ell}(\beta|\epsilon|/2)}
        {2\sinh(\beta\epsilon/2)} A(\epsilon),

    where :math:`i_\ell(x)` is the `modified spherical Bessel function of the first kind`_.

Norm of a solution :math:`A(\epsilon)` is defined as

    .. math::
        N = \int\limits_{-\infty}^\infty d\epsilon A(\epsilon) = \pi\chi(i\Omega_n=0).

The correlator of a real frequency is reconstructed according to

    .. math::
        \chi(\epsilon) = \frac{1}{\pi}\int\limits_{-\infty}^\infty
        d\epsilon' \frac{\epsilon' A(\epsilon')}{\epsilon' - \epsilon - i0}.

.. _bosonautocorr:

Autocorrelator of a boson-like operator
---------------------------------------

Enabled when `Som()` object is constructed with `kind = "BosonAutoCorr"`.

:math:`A(\epsilon)` is the spectral function to be found. It is defined differently from the fermionic
case, namely :math:`A(\epsilon) = \Im\chi(\epsilon)/\epsilon`, where

    .. math::
        \chi(\epsilon) = \int\limits_{-\infty}^\infty dt\ e^{i\epsilon t}\chi(t),\quad
        \chi(t) = i\theta(t)\langle[\hat O(t),\hat O(0)]\rangle.

Expressions in this section imply that :math:`A(-\epsilon) = A(\epsilon)`, and, therefore,
it is enough to find the spectral function on the positive energy half-axis only.
This condition is fulfilled by the dynamical susceptibilities of a form
:math:`\chi(\tau) = \langle\mathcal{T}\hat O(\tau)\hat O(0)\rangle`, where :math:`\hat O` is
a Hermitian operator.

- In imaginary time, :math:`\chi(\tau)`

    .. math::
        \chi(\tau) = \int\limits_0^\infty \frac{d\epsilon}{\pi}
        \frac{\epsilon (e^{-\tau\epsilon}+e^{-(\beta-\tau)\epsilon})}
        {1-e^{-\beta\epsilon}} A(\epsilon).

- At Matsubara frequencies, :math:`\chi(i\Omega_n)`

    .. math::
        \chi(i\Omega_n) = \int\limits_0^\infty
        d\epsilon \frac{1}{\pi}\frac{2\epsilon^2}{\Omega_n^2+\epsilon^2} A(\epsilon).

- In Legendre polynomial basis, :math:`\chi(\ell)`

    .. math::
        \chi(\ell) = \left\{
            \begin{array}{ll}
            \int\limits_0^\infty
            \frac{d\epsilon}{\pi}
            \frac{\beta\epsilon\sqrt{2\ell+1} i_{\ell}(\beta\epsilon/2)}
            {\sinh(\beta\epsilon/2)} A(\epsilon), &\ell\ \mathrm{ even},\\
            0, &\ell\ \mathrm{odd}.
        \end{array}\right.,

    where :math:`i_\ell(x)` is the `modified spherical Bessel function of the first kind`_.

Norm of a solution :math:`A(\epsilon)` is defined as

    .. math::
        N = \int\limits_0^\infty d\epsilon A(\epsilon) = \frac{\pi}{2}\chi(i\Omega_n=0).

The correlator of a real frequency is reconstructed according to

    .. math::
        \chi(\epsilon) =
        \frac{1}{\pi}\int\limits_{-\infty}^0
        d\epsilon' \frac{\epsilon' A(-\epsilon')}{\epsilon' - \epsilon - i0} +
        \frac{1}{\pi}\int\limits_0^\infty
        d\epsilon' \frac{\epsilon' A(\epsilon')}{\epsilon' - \epsilon - i0}.

.. _zerotemp:

Correlation function at zero temperature
----------------------------------------

:math:`A(\epsilon) = -(1/\pi)\Im G(\epsilon)` is the spectral function to be found.

Enabled when `Som()` object is constructed with `kind = "ZeroTemp"`.

- In imaginary time, :math:`G(\tau)`, :math:`\tau\in[0;\tau_{max}]`

    .. math::
        G(\tau) = -\int\limits_0^\infty
        d\epsilon\ e^{-\tau\epsilon} A(\epsilon).

- At Matsubara frequencies, :math:`G(i\omega_n)`, where :math:`\omega_n = (2n+1)\pi/\tau_{max}`
  for fermionic correlation functions and :math:`\omega_n = 2n\pi/\tau_{max}` for bosonic.

    .. math::
         G(i\omega_n) = \int\limits_0^\infty
         d\epsilon \frac{1}{i\omega_n-\epsilon} A(\epsilon).

- In Legendre polynomial basis, :math:`G(\ell)`

    .. math::
        G(\ell) = \int\limits_0^\infty
        d\epsilon \tau_{max}(-1)^{\ell+1}\sqrt{2\ell+1}
        i_{\ell}\left(\frac{\epsilon\tau_{max}}{2}\right)
        \exp\left(-\frac{\epsilon\tau_{max}}{2}\right) A(\epsilon),

    where :math:`i_\ell(x)` is the modified spherical Bessel function of the first kind.

Norm of a solution :math:`A(\epsilon)` is defined as

    .. math::
         N = \int\limits_0^\infty d\epsilon A(\epsilon).

The correlation function of a real frequency is reconstructed according to

    .. math::
         G(\epsilon) = -\int\limits_0^\infty
         d\epsilon' \frac{A(\epsilon')}{\epsilon' - \epsilon - i0}.
