.. _som_ref:

``som``: Main Python module of SOM
----------------------------------

.. currentmodule:: som

``som.Som``
~~~~~~~~~~~

.. autoclass:: Som

    .. automethod:: __init__
    .. automethod:: accumulate(**kwargs)
    .. autoattribute:: last_accumulate_parameters
    .. autoattribute:: accumulate_status
    .. automethod:: adjust_f(**kwrags)
    .. automethod:: compute_final_solution(good_chi_rel: float = 2.0, good_chi_abs: float = cmath.inf, verbosity: int = 0)
    .. automethod:: compute_final_solution_cc(**kwargs)
    .. automethod:: run(**kwargs)
    .. automethod:: clear()
    .. autoattribute:: observable_kind
    .. autoattribute:: dim
    .. automethod:: particular_solutions(i: int)
    .. autoattribute:: objf_min
    .. automethod:: solution(i: int)
    .. autoattribute:: solutions
    .. automethod:: objf(i: int)
    .. autoattribute:: objf_list
    .. automethod:: histogram(i: int)
    .. autoattribute:: histograms

Free functions
~~~~~~~~~~~~~~

.. autofunction:: fill_refreq
.. autofunction:: compute_tail
.. autofunction:: reconstruct

.. autofunction:: estimate_boson_corr_spectrum_norms
.. autofunction:: count_good_solutions

Auxiliary types
~~~~~~~~~~~~~~~

.. autoclass:: Rectangle

    .. autoattribute:: center
    .. autoattribute:: width
    .. autoattribute:: height
    .. autoattribute:: norm
    .. automethod:: __call__(epsilon: float)

.. autoclass:: Configuration

    .. automethod:: __len__
    .. automethod:: __getitem__
    .. automethod:: __iter__
    .. automethod:: __call__(epsilon: float)
