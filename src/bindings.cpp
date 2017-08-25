#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <blitz/array.h>
#include <limits.h>

#include "cloud_identification.h"
#include "common.h"

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


py::array_t<indexint> number_objects(py::array_t<double> scalar_field, py::array_t<bool> mask)
{
  const size_t ndim = 3;

  py::buffer_info info_scalar_field = scalar_field.request();
  py::buffer_info info_mask = mask.request();

  if (info_scalar_field.ndim != ndim || info_mask.ndim != ndim) {
    throw std::runtime_error("Inputs should be 3D");
  }

  for (int n=0; n < info_scalar_field.ndim; n++) {
    if (info_scalar_field.shape[n] != info_mask.shape[n]) {
      throw std::runtime_error("Input shapes must be equal");
    }
  }

  py::array_t<indexint> result = py::array(py::buffer_info(
    nullptr,            /* Pointer to data (nullptr -> ask NumPy to allocate!) */
    sizeof(indexint),     /* Size of one item */
    info_scalar_field.format, /* Buffer format */
    info_scalar_field.ndim,          /* How many dimensions? */
    info_scalar_field.shape,  /* Number of elements for each dimension */
    info_scalar_field.strides  /* Strides for each dimension */
  ));


  blitz::Array<double,ndim> scalar_field_blitz = py_array_to_blitz<double,ndim>(scalar_field);
  blitz::Array<bool,ndim> mask_blitz = py_array_to_blitz<bool,ndim>(mask);
  blitz::Array<indexint,ndim> result_blitz = py_array_to_blitz<indexint,ndim>(result);

  double s_min = min(scalar_field_blitz);
  double s_max = max(scalar_field_blitz);
  double s_range = s_max - s_min;
  double ds = s_range/double(SHRT_MAX - SHRT_MIN);
  blitz::TinyVector<int,ndim> shape = scalar_field_blitz.shape();

  blitz::Array<short,ndim> rescaled_scalar_field_blitz = blitz::Array<short,ndim>(scalar_field_blitz.shape());
  for (int i=0; i<shape[0]; i++) {
    for (int j=0; j<shape[1]; j++) {
      for (int k=0; k<shape[2]; k++) {
        rescaled_scalar_field_blitz(i,j,k) = SHRT_MIN + (short)((scalar_field_blitz(i,j,k) - s_min)/ds);
      }
    }
  }

  find_objects(rescaled_scalar_field_blitz, mask_blitz, result_blitz);

  return result;
}

PYBIND11_PLUGIN(cloud_identification)
{
    py::module m("cloud_identification");
    m.def("number_objects", &number_objects, "Identify individual cloud objects in regions defined by mask",
          py::arg("scalar_field"), py::arg("mask"));
    return m.ptr();
}
