// DO NOT EDIT
// Generated automatically using libclang using the command :
// c++2py.py ../c++/som_core.hpp --only_converters --compiler_options=-DCACHE_SIZE=0x10000 -p -m som -o som --appname triqs_som --moduledoc "The Stochastic Optimization Method"


// --- C++ Python converter for run_parameters_t

namespace triqs { namespace py_tools {

template <> struct py_converter<run_parameters_t> {
 static PyObject *c2py(run_parameters_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "energy_window"      , convert_to_python(x.energy_window));
  PyDict_SetItemString( d, "max_time"           , convert_to_python(x.max_time));
  PyDict_SetItemString( d, "verbosity"          , convert_to_python(x.verbosity));
  PyDict_SetItemString( d, "t"                  , convert_to_python(x.t));
  PyDict_SetItemString( d, "f"                  , convert_to_python(x.f));
  PyDict_SetItemString( d, "adjust_f"           , convert_to_python(x.adjust_f));
  PyDict_SetItemString( d, "l"                  , convert_to_python(x.l));
  PyDict_SetItemString( d, "adjust_l"           , convert_to_python(x.adjust_l));
  PyDict_SetItemString( d, "make_histograms"    , convert_to_python(x.make_histograms));
  PyDict_SetItemString( d, "random_seed"        , convert_to_python(x.random_seed));
  PyDict_SetItemString( d, "random_name"        , convert_to_python(x.random_name));
  PyDict_SetItemString( d, "max_rects"          , convert_to_python(x.max_rects));
  PyDict_SetItemString( d, "min_rect_width"     , convert_to_python(x.min_rect_width));
  PyDict_SetItemString( d, "min_rect_weight"    , convert_to_python(x.min_rect_weight));
  PyDict_SetItemString( d, "distrib_d_max"      , convert_to_python(x.distrib_d_max));
  PyDict_SetItemString( d, "gamma"              , convert_to_python(x.gamma));
  PyDict_SetItemString( d, "adjust_f_range"     , convert_to_python(x.adjust_f_range));
  PyDict_SetItemString( d, "adjust_f_l"         , convert_to_python(x.adjust_f_l));
  PyDict_SetItemString( d, "adjust_f_kappa"     , convert_to_python(x.adjust_f_kappa));
  PyDict_SetItemString( d, "adjust_l_range"     , convert_to_python(x.adjust_l_range));
  PyDict_SetItemString( d, "adjust_l_good_d"    , convert_to_python(x.adjust_l_good_d));
  PyDict_SetItemString( d, "adjust_l_verygood_d", convert_to_python(x.adjust_l_verygood_d));
  PyDict_SetItemString( d, "adjust_l_ratio"     , convert_to_python(x.adjust_l_ratio));
  PyDict_SetItemString( d, "hist_max"           , convert_to_python(x.hist_max));
  PyDict_SetItemString( d, "hist_n_bins"        , convert_to_python(x.hist_n_bins));
  return d;
 }

 template <typename T, typename U> static void _get_optional(PyObject *dic, const char *name, T &r, U const &init_default) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = init_default;
 }

 template <typename T> static void _get_optional(PyObject *dic, const char *name, T &r) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = T{};
 }

 static run_parameters_t py2c(PyObject *dic) {
  run_parameters_t res;
  res.energy_window = convert_from_python<std::pair<double, double>>(PyDict_GetItemString(dic, "energy_window"));
  _get_optional(dic, "max_time"           , res.max_time              ,-1);
  _get_optional(dic, "verbosity"          , res.verbosity             ,((triqs::mpi::communicator().rank()==0)?2:0));
  _get_optional(dic, "t"                  , res.t                     ,1000);
  _get_optional(dic, "f"                  , res.f                     ,100);
  _get_optional(dic, "adjust_f"           , res.adjust_f              ,true);
  _get_optional(dic, "l"                  , res.l                     ,500);
  _get_optional(dic, "adjust_l"           , res.adjust_l              ,false);
  _get_optional(dic, "make_histograms"    , res.make_histograms       ,false);
  _get_optional(dic, "random_seed"        , res.random_seed           ,34788+928374*triqs::mpi::communicator().rank());
  _get_optional(dic, "random_name"        , res.random_name           ,"");
  _get_optional(dic, "max_rects"          , res.max_rects             ,60);
  _get_optional(dic, "min_rect_width"     , res.min_rect_width        ,1e-3);
  _get_optional(dic, "min_rect_weight"    , res.min_rect_weight       ,1e-3);
  _get_optional(dic, "distrib_d_max"      , res.distrib_d_max         ,2);
  _get_optional(dic, "gamma"              , res.gamma                 ,2);
  _get_optional(dic, "adjust_f_range"     , res.adjust_f_range        ,std::pair<int,int>{100,5000});
  _get_optional(dic, "adjust_f_l"         , res.adjust_f_l            ,20);
  _get_optional(dic, "adjust_f_kappa"     , res.adjust_f_kappa        ,0.25);
  _get_optional(dic, "adjust_l_range"     , res.adjust_l_range        ,std::pair<int,int>{100,2000});
  _get_optional(dic, "adjust_l_good_d"    , res.adjust_l_good_d       ,2.0);
  _get_optional(dic, "adjust_l_verygood_d", res.adjust_l_verygood_d   ,4.0/3.0);
  _get_optional(dic, "adjust_l_ratio"     , res.adjust_l_ratio        ,0.95);
  _get_optional(dic, "hist_max"           , res.hist_max              ,2.0);
  _get_optional(dic, "hist_n_bins"        , res.hist_n_bins           ,100);
  return res;
 }

 template <typename T>
 static void _check(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!convertible_from_python<T>(PyDict_GetItemString(dic, name), false))
   fs << "\n" << ++err << " The parameter " << name << " does not have the right type : expecting " << tname
      << " in C++, but got '" << PyDict_GetItemString(dic, name)->ob_type->tp_name << "' in Python.";
 }

 template <typename T>
 static void _check_mandatory(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!PyDict_Contains(dic, pyref::string(name)))
   fs << "\n" << ++err << " Mandatory parameter " << name << " is missing.";
  else _check<T>(dic,fs,err,name,tname);
 }

 template <typename T>
 static void _check_optional(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (PyDict_Contains(dic, pyref::string(name))) _check<T>(dic, fs, err, name, tname);
 }

 static bool is_convertible(PyObject *dic, bool raise_exception) {
  if (!PyDict_Check(dic)) {
   if (raise_exception) { PyErr_SetString(PyExc_TypeError, "Not a python dict");}
   return false;
  }
  std::stringstream fs, fs2; int err=0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
  std::vector<std::string> ks, all_keys = {"energy_window","max_time","verbosity","t","f","adjust_f","l","adjust_l","make_histograms","random_seed","random_name","max_rects","min_rect_width","min_rect_weight","distrib_d_max","gamma","adjust_f_range","adjust_f_l","adjust_f_kappa","adjust_l_range","adjust_l_good_d","adjust_l_verygood_d","adjust_l_ratio","hist_max","hist_n_bins"};
  pyref keys = PyDict_Keys(dic);
  if (!convertible_from_python<std::vector<std::string>>(keys, true)) {
   fs << "\nThe dict keys are not strings";
   goto _error;
  }
  ks = convert_from_python<std::vector<std::string>>(keys);
  for (auto & k : ks)
   if (std::find(all_keys.begin(), all_keys.end(), k) == all_keys.end())
    fs << "\n"<< ++err << " The parameter '" << k << "' is not recognized.";
#endif

  _check_mandatory<std::pair<double, double>>(dic, fs, err, "energy_window"      , "std::pair<double, double>");
  _check_optional <int                      >(dic, fs, err, "max_time"           , "int");
  _check_optional <int                      >(dic, fs, err, "verbosity"          , "int");
  _check_optional <int                      >(dic, fs, err, "t"                  , "int");
  _check_optional <int                      >(dic, fs, err, "f"                  , "int");
  _check_optional <bool                     >(dic, fs, err, "adjust_f"           , "bool");
  _check_optional <int                      >(dic, fs, err, "l"                  , "int");
  _check_optional <bool                     >(dic, fs, err, "adjust_l"           , "bool");
  _check_optional <bool                     >(dic, fs, err, "make_histograms"    , "bool");
  _check_optional <int                      >(dic, fs, err, "random_seed"        , "int");
  _check_optional <std::string              >(dic, fs, err, "random_name"        , "std::string");
  _check_optional <int                      >(dic, fs, err, "max_rects"          , "int");
  _check_optional <double                   >(dic, fs, err, "min_rect_width"     , "double");
  _check_optional <double                   >(dic, fs, err, "min_rect_weight"    , "double");
  _check_optional <double                   >(dic, fs, err, "distrib_d_max"      , "double");
  _check_optional <double                   >(dic, fs, err, "gamma"              , "double");
  _check_optional <std::pair<int, int>      >(dic, fs, err, "adjust_f_range"     , "std::pair<int, int>");
  _check_optional <int                      >(dic, fs, err, "adjust_f_l"         , "int");
  _check_optional <double                   >(dic, fs, err, "adjust_f_kappa"     , "double");
  _check_optional <std::pair<int, int>      >(dic, fs, err, "adjust_l_range"     , "std::pair<int, int>");
  _check_optional <double                   >(dic, fs, err, "adjust_l_good_d"    , "double");
  _check_optional <double                   >(dic, fs, err, "adjust_l_verygood_d", "double");
  _check_optional <double                   >(dic, fs, err, "adjust_l_ratio"     , "double");
  _check_optional <double                   >(dic, fs, err, "hist_max"           , "double");
  _check_optional <int                      >(dic, fs, err, "hist_n_bins"        , "int");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class run_parameters_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}
