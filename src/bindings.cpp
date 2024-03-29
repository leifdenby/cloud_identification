#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <blitz/array.h>
#include <limits.h>
#include <cmath>

#include "cloud_identification.h"
#include "minkowski.h"
#include "common.h"
#include "filters.h"

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
py::array_t<T> blitz_array_to_py(blitz::Array<T,N> &a) {
  std::vector<py::size_t> shape = std::vector<py::size_t>(N);
  std::vector<py::size_t> strides = std::vector<py::size_t>(N);

  for (int i = 0; i < N; i++)
  {
    shape[i] = a.shape()[i];
    strides[i] = a.stride()[i]*sizeof(T);
  }

  py::array_t<T> result = py::array(py::buffer_info(
    a.data(),
    sizeof(T),
    py::format_descriptor<T>::format(),
    a.dimensions(),
    shape,
    strides
  ));

  return result;
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

  blitz::TinyVector<int,ndim> shape = scalar_field_blitz.shape();
  blitz::Array<short,ndim> rescaled_scalar_field_blitz = blitz::Array<short,ndim>(scalar_field_blitz.shape());

  scale_field(scalar_field_blitz, rescaled_scalar_field_blitz);

  find_objects(rescaled_scalar_field_blitz, mask_blitz, result_blitz);

  return result;
}


py::array_t<short> apply_scaling(py::array_t<double> scalar_field)
{
  const size_t ndim = 3;

  py::buffer_info info_scalar_field = scalar_field.request();

  py::array_t<short> result = py::array(py::buffer_info(
    nullptr,            /* Pointer to data (nullptr -> ask NumPy to allocate!) */
    sizeof(indexint),     /* Size of one item */
    info_scalar_field.format, /* Buffer format */
    info_scalar_field.ndim,          /* How many dimensions? */
    info_scalar_field.shape,  /* Number of elements for each dimension */
    info_scalar_field.strides  /* Strides for each dimension */
  ));

  blitz::Array<double,ndim> scalar_field_blitz = py_array_to_blitz<double,ndim>(scalar_field);
  blitz::Array<short,ndim> result_blitz = py_array_to_blitz<short,ndim>(result);

  scale_field(scalar_field_blitz, result_blitz);

  return result;
}



void remove_intersecting(py::array_t<indexint> labels, py::array_t<bool> mask)
{
  const size_t ndim = 3;

  py::buffer_info info_labels = labels.request();
  py::buffer_info info_mask = mask.request();

  if (info_labels.ndim != ndim || info_mask.ndim != ndim) {
    throw std::runtime_error("Inputs should be 3D");
  }

  for (int n=0; n < info_labels.ndim; n++) {
    if (info_labels.shape[n] != info_mask.shape[n]) {
      throw std::runtime_error("Input shapes must be equal");
    }
  }

  blitz::Array<indexint,ndim> labels_blitz = py_array_to_blitz<indexint,ndim>(labels);
  blitz::Array<bool,ndim> mask_blitz = py_array_to_blitz<bool,ndim>(mask);

  filters::remove_intersecting(labels_blitz, mask_blitz);
}


void filter_labels(py::array_t<indexint> labels, py::array_t<indexint> idxs_keep)
{
  const size_t ndim = 3;

  py::buffer_info info_labels = labels.request();
  py::buffer_info info_idxs_keep = idxs_keep.request();

  if (info_labels.ndim != ndim) {
    throw std::runtime_error("Labels should be 3D");
  }

  if (info_idxs_keep.ndim != 1) {
    throw std::runtime_error("Labels to keep should be 1D");
  }

  blitz::Array<indexint,ndim> labels_blitz = py_array_to_blitz<indexint,ndim>(labels);
  blitz::Array<indexint,1> idxs_keep_blitz = py_array_to_blitz<indexint,1>(idxs_keep);

  filters::filter_labels(labels_blitz, idxs_keep_blitz);
}


void remove_intersecting_xy(py::array_t<indexint> labels, py::array_t<bool> mask)
{
  const size_t ndim = 3;

  py::buffer_info info_labels = labels.request();
  py::buffer_info info_mask = mask.request();

  if (info_labels.ndim != ndim || info_mask.ndim != ndim-1) {
    throw std::runtime_error("Inputs should be 3D for object labels and 2D for mask");
  }

  for (int n=0; n < info_labels.ndim-1; n++) {
    if (info_labels.shape[n] != info_mask.shape[n]) {
      throw std::runtime_error("Input shapes must be equal");
    }
  }

  blitz::Array<indexint,ndim> labels_blitz = py_array_to_blitz<indexint,ndim>(labels);
  blitz::Array<bool,ndim-1> mask_blitz = py_array_to_blitz<bool,ndim-1>(mask);

  filters::remove_intersecting_xy(labels_blitz, mask_blitz);
}


py::array_t<int> N0(py::array_t<int> labels)
{
  const size_t ndim = 3;
  py::buffer_info info_labels = labels.request();
  if (info_labels.ndim != ndim) {
    throw std::runtime_error("Input should be 3D");
  }

  blitz::Array<int,ndim> labels_blitz = py_array_to_blitz<int,ndim>(labels);
  blitz::Array<int,1> n_vertices = minkowski::N0(labels_blitz);

  return blitz_array_to_py<int,1>(n_vertices);
}

py::array_t<int> N1(py::array_t<int> labels)
{
  const size_t ndim = 3;
  py::buffer_info info_labels = labels.request();
  if (info_labels.ndim != ndim) {
    throw std::runtime_error("Input should be 3D");
  }

  blitz::Array<int,ndim> labels_blitz = py_array_to_blitz<int,ndim>(labels);
  blitz::Array<int,1> n_edges = minkowski::N1(labels_blitz);

  return blitz_array_to_py<int,1>(n_edges);
}

py::array_t<int> N2(py::array_t<int> labels)
{
  const size_t ndim = 3;
  py::buffer_info info_labels = labels.request();
  if (info_labels.ndim != ndim) {
    throw std::runtime_error("Input should be 3D");
  }

  blitz::Array<int,ndim> labels_blitz = py_array_to_blitz<int,ndim>(labels);
  blitz::Array<int,1> n_faces = minkowski::N2(labels_blitz);

  return blitz_array_to_py<int,1>(n_faces);
}

py::array_t<int> N3(py::array_t<int> labels)
{
  const size_t ndim = 3;
  py::buffer_info info_labels = labels.request();
  if (info_labels.ndim != ndim) {
    throw std::runtime_error("Input should be 3D");
  }

  blitz::Array<int,ndim> labels_blitz = py_array_to_blitz<int,ndim>(labels);
  blitz::Array<int,1> n_cells = minkowski::N3(labels_blitz);

  return blitz_array_to_py<int,1>(n_cells);
}

py::array_t<float> topological_scales(py::array_t<int> labels, float dx)
{
  const size_t ndim = 3;
  py::buffer_info info_labels = labels.request();
  if (info_labels.ndim != ndim) {
    throw std::runtime_error("Input should be 3D");
  }

  blitz::Array<int,ndim> labels_blitz = py_array_to_blitz<int,ndim>(labels);
  blitz::Array<float,2> scales = minkowski::topological_scales(labels_blitz, dx);

  return blitz_array_to_py<float,2>(scales);
}


py::array_t<float> filamentarity(py::array_t<float> mf)
{
  const size_t ndim = 2;
  py::buffer_info info_labels = mf.request();
  if (info_labels.ndim != ndim) {
    throw std::runtime_error("Input should be 2D");
  }

  blitz::Array<float,ndim> mf_blitz = py_array_to_blitz<float,ndim>(mf);
  blitz::Array<float,1> filamentarity = minkowski::filamentarity(mf_blitz);

  return blitz_array_to_py<float,1>(filamentarity);
}


py::array_t<float> planarity(py::array_t<float> mf)
{
  const size_t ndim = 2;
  py::buffer_info info_labels = mf.request();
  if (info_labels.ndim != ndim) {
    throw std::runtime_error("Input should be 2D");
  }

  blitz::Array<float,ndim> mf_blitz = py_array_to_blitz<float,ndim>(mf);
  blitz::Array<float,1> planarity = minkowski::planarity(mf_blitz);

  return blitz_array_to_py<float,1>(planarity);
}

PYBIND11_PLUGIN(cloud_identification)
{
    py::module m("cloud_identification");
    m.def("number_objects", &number_objects, "Identify individual cloud objects in regions defined by mask",
          py::arg("scalar_field"), py::arg("mask"));
    m.def("apply_scaling", &apply_scaling, "Apply scaling which is used internally when number objects",
          py::arg("scalar_field"));
    m.def("remove_intersecting", &remove_intersecting, "Remove all labelled objects which intersect with mask",
          py::arg("labels"), py::arg("mask"));
    m.def("remove_intersecting_xy", &remove_intersecting, "Remove all labelled objects which intersect with 2D xy mask",
          py::arg("labels"), py::arg("mask"));
    m.def("filter_labels", &filter_labels, "Remove labelled objects which aren't in list of indexes",
          py::arg("labels"), py::arg("idxs_keep"));
    m.def("N0", &N0, "Find number of vertices for each labelled object");
    m.def("N1", &N1, "Find number of edges for each labelled object");
    m.def("N2", &N2, "Find number of faces for each labelled object");
    m.def("N3", &N3, "Find number of cubes (i.e the volume) for each labelled object");
    m.def("topological_scales", &topological_scales, "Compute characteristic topological scales (in terms of characteristic thickness, width, length and genus) using Minkowski functionals", 
          py::arg("labels"), py::arg("dx"));
    m.def("planarity", &planarity, "Compute planarity of objects from Minkowski functionals", py::arg("mf"));
    m.def("filamentarity", &filamentarity, "Compute filamentarity of objects from Minkowski functionals", py::arg("mf"));
    return m.ptr();
}
