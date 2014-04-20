from libcpp.string cimport string
from arrays cimport *
from pytriqs.gf.local.gf cimport *
import numpy
from pytriqs.gf.local import *
import pytriqs.utility.mpi as mpi
from pytriqs.parameters.parameters cimport *

cdef extern from "c++/som.hpp" :

    # Declare the C++ class
    cdef cppclass my_class_c "my_class":
      
      my_class_c(parameters) except +

      #input containers
      double beta()
      gf_block_imtime & g1_view()
      gf_refreq & g2_view()
      gf_imtime & g3_view()
      array_view Arr_view()
      array_view ar_view()
      matrix_view[double] m1_view()

      parameter_defaults param_defaults() 
  
cdef class MyClass:

    cdef my_class_c * _c

    def __init__(self, **kw):
        assert ('parameters' in kw and len(kw)==1) or ('parameters' not in kw), 'wrong arguments'
        cdef object param
        param=kw['parameters'] if 'parameters' in kw else Parameters().update(kw)
        self._c = new my_class_c((<Parameters?>param)._c)

    def __dealloc__(self):
        del self._c

    property beta:
        """ Doc """
        def __get__(self): 
            return self._c.beta()

    property g1:
        """ Doc of g1 ?"""
        def __get__(self): 
            return make_BlockGfImTime(self._c.g1_view(), [['0'],['0']], "g1(tau)")

    property g2:
        """ Doc of g2 ?"""
        def __get__(self): 
            return make_GfReFreq(self._c.g2_view(), [['0'],['0']], "g2(w)")

    property g3:
        """ Doc of g3 ?"""
        def __get__(self): 
            return make_GfImTime(self._c.g3_view(), [['0'],['0']], "g3(tau)")

    property Arr:
        """ Doc of Arr ?"""
        def __get__(self): 
            return (self._c.Arr_view().to_python(), "Arr")

    property ar:
        """ Doc of ar ?"""
        def __get__(self): 
            return (self._c.ar_view().to_python(), "ar")

    property m:
        """ Doc of m ?"""
        def __get__(self): 
            return self._c.m1_view().to_python()

