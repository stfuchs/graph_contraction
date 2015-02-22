/**
 * @file   graphcontraction.cpp
 * @author Steffen Fuchs <steffen@steffen-ubuntu>
 * @date   Sun Feb 22 14:18:37 2015
 * 
 * @brief  
 * 
 * 
 */


#include "gc/core.hpp"
#include "gc/policies.hpp"
#include <boost/python.hpp>

#define BIND_GC(type, name)\
{\
class_<type>(name, init<float,int>())\
  .def("fit", &type::fit)\
  .def("set_grid_adjacency", &type::set_grid_adjacency)\
  .def("get_labels", &type::get_labels)\
  .def("get_representer", &type::get_representer);\
}


using namespace boost::python;
typedef GC::GraphContraction<double, GC::VariancePolicy> VGC_double;
typedef GC::GraphContraction<float, GC::VariancePolicy> VGC_float;


BOOST_PYTHON_MODULE(graphcontraction)
{
  BIND_GC(VGC_double, "VarianceGC");
  BIND_GC(VGC_float , "VarianceGC");
}
