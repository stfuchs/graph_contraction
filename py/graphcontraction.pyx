from libcpp.vector cimport vector
from libcpp.pair cimport pair
#cimport libeigenmap as lem
from libeigenmap cimport MapXi, MapXd, MapX3f, MapX3d
cimport numpy as np
import numpy as np
import ctypes

cdef extern from "../include/gc/quad.hpp":
    cdef cppclass QGC_VAR_F3 "GC::QuadHierarchicalContraction<float,3>":
        QGC_VAR_F3(float,float,int) except+
        void init_data(int,int,MapX3f&)
        void fit()
        void get_labels(MapXi&)
        void get_representer(MapX3f&)

cdef class QGC_F3:
    """Implements GC::QuadHierarichalContraction<float,3> from quad.hpp"""
    cdef QGC_VAR_F3* ptr
    cdef int ndata
    def __cinit__(self, float ratio, float max_var, int level):
        self.ptr = new QGC_VAR_F3(ratio, max_var, level)
    def __dealloc__(self):
        del self.ptr
    cpdef fit(self, int height, int width, object[float, ndim=2] data):
        if data.shape[1] != 3:
            raise ValueError(
                "This is the 3D implementation! Your dimension is %s"%data.shape[1])
        cdef MapX3f* m = new MapX3f(&data[0,0],data.shape[0])
        self.ptr.init_data(height,width,m[0])
        self.ptr.fit()
        self.ndata = data.shape[0]
        del m
    cpdef get_labels(self):
        cdef np.ndarray[int,ndim=1] res = np.zeros(self.ndata,dtype=ctypes.c_int)
        cdef MapXi* m = new MapXi(&res[0],self.ndata)
        self.ptr.get_labels(m[0])
        del m
        return res
    cpdef get_representer(self):
        cdef np.ndarray[float,ndim=2] res = np.zeros([self.ndata,3],dtype=ctypes.c_float)
        cdef MapX3f* m = new MapX3f(&res[0,0],self.ndata)
        self.ptr.get_representer(m[0])
        del m
        return res

