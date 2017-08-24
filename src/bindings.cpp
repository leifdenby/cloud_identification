#include "math.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <blitz/array.h>

namespace py = pybind11;

template<class T, int N>
blitz::Array<T,N> py_array_to_blitz(py::array_t<T> &a) {
  blitz::TinyVector<int,N> shape(0);
  blitz::TinyVector<int,N> strides(0);

  py::buffer_info info_a = a.request();

  for (int i = 0; i < N; i++)
  {
    shape[i] = info_a.shape[i];
    strides[i] = info_a.strides[i]/sizeof(T);
  }

  return blitz::Array<T,N>((T*) info_a.ptr, shape, strides, blitz::neverDeleteData);
}

template<class T, int N>
void blitz_add(
  const blitz::Array<T,N> &a,
  const blitz::Array<T,N> &b,
  blitz::Array<T,N> &c
) {

  if (N == 1) {
    for (int i=0; i < a.size(); i++) {
      c(i) = a(i) + b(i);
    }
  }
  else if (N == 2) {
    for (int i=0; i < a.shape()[0]; i++) {
      for (int j=0; j < a.shape()[1]; j++) {
        c(i,j) = a(i,j) + b(i,j);
      }
    }
  }
}


py::array_t<double> add_arrays(py::array_t<double> a, py::array_t<double> b)
{
  const size_t ndim = 2;

  py::buffer_info info_a = a.request();
  py::buffer_info info_b = b.request();

  if (info_a.ndim != ndim || info_b.ndim != ndim) {
    throw std::runtime_error("Inputs should be 2D");
  }

  for (int n=0; n < info_a.ndim; n++) {
    if (info_a.shape[n] != info_b.shape[n]) {
      throw std::runtime_error("Input shapes must be equal");
    }
  }

  py::array_t<double> result = py::array(py::buffer_info(
    nullptr,            /* Pointer to data (nullptr -> ask NumPy to allocate!) */
    sizeof(double),     /* Size of one item */
    info_a.format, /* Buffer format */
    info_a.ndim,          /* How many dimensions? */
    info_a.shape,  /* Number of elements for each dimension */
    info_a.strides  /* Strides for each dimension */
  ));


  blitz::Array<double,ndim> a_blitz = py_array_to_blitz<double,ndim>(a);
  blitz::Array<double,ndim> b_blitz = py_array_to_blitz<double,ndim>(b);
  blitz::Array<double,ndim> c_blitz = py_array_to_blitz<double,ndim>(result);

  blitz_add<double,ndim>(a_blitz, b_blitz, c_blitz);

  return result;
}    

PYBIND11_PLUGIN(python_cpp_example)
{
    py::module m("python_cpp_example");
    m.def("add_arrays", &add_arrays, "Adding two numpy arrays");

    m.def("add", &add);
    m.def("subtract", &subtract);
    return m.ptr();
}
