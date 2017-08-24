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

  for (int i=0; i < a.size(); i++) {
    c(i) = a(i) + b(i);
  }
}


py::array_t<double> add_arrays(py::array_t<double> a, py::array_t<double> b)
{
  py::buffer_info info_a = a.request();
  py::buffer_info info_b = b.request();

  if ((info_a.ndim != 1) || (info_b.ndim != 1))
    throw std::runtime_error("Number of dimensions must be one");

  if (info_a.shape[0] != info_b.shape[0])
    throw std::runtime_error("Input shapes must be equal");

  //std::vector<double> ret(info_a.shape[0]);
  //for (unsigned int idx = 0; idx < info_a.shape[0]; idx++)
    //ret[idx] = ((double*)info_a.ptr)[idx] + ((double*)info_b.ptr)[idx];
    //

  py::array_t<double> result = py::array_t<double>(info_a.size);

  std::vector<size_t> strides = {sizeof(double)};
  std::vector<size_t> shape = {info_a.shape[0]};
  const size_t ndim = 1;

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
