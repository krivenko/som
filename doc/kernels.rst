.. _kernels:

Supported integral kernels
==========================

.. note::

    :math:`\beta` denotes the inverse temperature in all equations below.
    :math:`A(\epsilon)` is the spectral function to be found.

Fermionic thermal Green's function
----------------------------------

Enabled when `Som()` object is constucted with `kind = "FermionGf"`.

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

  where :math:`i_\ell(x)` is the modified spherical Bessel function of the first kind.
