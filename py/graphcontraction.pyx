from libcpp.vector cimport vector
from libcpp.pair cimport pair
#cimport libeigenmap as lem
from libeigenmap cimport MapXi, MapXd, MapX3f, MapX3d
cimport numpy as np
import numpy as np
import ctypes

cdef extern from "../include/gc/core.hpp":
    cdef cppclass GC_VAR_D1C "GC::GraphContraction<double,1,GC::Variance>":
        GC_VAR_D1C(double) except+
        void init_grid_adjacency(int,int)
        void init_adjacency(vector[pair[int,int]]&, int)
        void init_data(MapXd&)
        void fit()
        void get_labels(MapXi&)
        void get_representer(MapXd&)
        int data_size()
    cdef cppclass GC_VAR_F3C "GC::GraphContraction<float,3,GC::Variance>":
        GC_VAR_F3C(float) except+
        void init_grid_adjacency(int,int)
        void init_adjacency(vector[pair[int,int]]&, int)
        void init_data(MapX3f&)
        void fit()
        void get_labels(MapXi&)
        void get_representer(MapX3f&)
        int data_size()
    cdef cppclass GC_VAR_D3C "GC::GraphContraction<double,3,GC::Variance>":
        GC_VAR_D3C(double) except+
        void init_grid_adjacency(int,int)
        void init_adjacency(vector[pair[int,int]]&, int)
        void init_data(MapX3d&)
        void fit()
        void get_labels(MapXi&)
        void get_representer(MapX3d&)
        int data_size()
cdef extern from "../include/gc/quad.hpp":
    cdef cppclass QGC_VAR_F3 "GC::QuadHierarchicalContraction<float,3>":
        QGC_VAR_F3(float,int) except+
        void init_data(int,int,MapX3f&)
        void fit()
        void get_labels(MapXi&)
        void get_representer(MapX3f&)

cdef class GC_D1:
    """Implements GC::Graphcontraction<double,1,GC::Variance> from core.hpp"""
    cdef GC_VAR_D1C* ptr
    cdef int ndata
    def __cinit__(self, double max_var):
        self.ptr = new GC_VAR_D1C(max_var)
    def __dealloc__(self):
        del self.ptr
    def init_grid_adjacency(self, int height, int width):
        """ set graph adjacency as grid (height, width) """
        self.ptr.init_grid_adjacency(height, width)
    def init_adjacency(self, edges, ndata):
        """ set graph adjacency manually
        edges: list of pairs of vertex indices:
               [(v_idx1,v_idx2),(idx...,idx...),...]
        ndata: number of data points
        """
        self.ptr.init_adjacency(edges, ndata)
    cpdef fit(self, object[np.double_t, ndim=1] data):
        if self.ptr.data_size() != data.size:
            raise ValueError(
                "data.size does not match adjacency. Expected %s"%self.ptr.data_size())
        cdef MapXd* m = new MapXd(&data[0],data.size)
        self.ptr.init_data(m[0])
        self.ptr.fit()
        self.ndata = data.size
        del m
    cpdef get_labels(self):
        cdef np.ndarray[int,ndim=1] res = np.zeros(self.ndata,dtype=ctypes.c_int)
        cdef MapXi* m = new MapXi(&res[0],self.ndata)
        self.ptr.get_labels(m[0])
        del m
        return res
    cpdef get_representer(self):
        cdef np.ndarray[np.double_t,ndim=1] res = np.zeros(self.ndata,dtype=np.double)
        cdef MapXd* m = new MapXd(&res[0],self.ndata)
        self.ptr.get_representer(m[0])
        del m
        return res


cdef class GC_F3:
    """Implements GC::Graphcontraction<float,3,GC::Variance> from core.hpp"""
    cdef GC_VAR_F3C* ptr
    cdef int ndata
    def __cinit__(self, float max_var):
        self.ptr = new GC_VAR_F3C(max_var)
    def __dealloc__(self):
        del self.ptr
    def init_grid_adjacency(self, int height, int width):
        """ set graph adjacency as grid (height, width) """
        self.ptr.init_grid_adjacency(height, width)
    def init_adjacency(self, edges, ndata):
        """ set graph adjacency manually
        edges: list of pairs of vertex indices:
               [(v_idx1,v_idx2),(idx...,idx...),...]
        ndata: number of data points
        """
        self.ptr.init_adjacency(edges, ndata)
    cpdef fit(self, object[float, ndim=2] data):
        if data.shape[1] != 3:
            raise ValueError(
                "This is the 3D implementation! Your dimension is %s"%data.shape[1])
        if self.ptr.data_size() != data.shape[0]:
            raise ValueError(
                "data.size does not match adjacency. Expected %s"%self.ptr.data_size())
        cdef MapX3f* m = new MapX3f(&data[0,0],data.shape[0])
        self.ptr.init_data(m[0])
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


cdef class GC_D3:
    """Implements GC::Graphcontraction<double,3,GC::Variance> from core.hpp"""
    cdef GC_VAR_D3C* ptr
    cdef int ndata
    def __cinit__(self, double max_var):
        self.ptr = new GC_VAR_D3C(max_var)
    def __dealloc__(self):
        del self.ptr
    def init_grid_adjacency(self, int height, int width):
        """ set graph adjacency as grid (height, width) """
        self.ptr.init_grid_adjacency(height, width)
    def init_adjacency(self, edges, ndata):
        """ set graph adjacency manually
        edges: list of pairs of vertex indices:
               [(v_idx1,v_idx2),(idx...,idx...),...]
        ndata: number of data points
        """
        self.ptr.init_adjacency(edges, ndata)
    cpdef fit(self, object[np.double_t, ndim=2] data):
        if data.shape[1] != 3:
            raise ValueError(
                "This is the 3D implementation! Your dimension is %s"%data.shape[1])
        if self.ptr.data_size() != data.shape[0]:
            raise ValueError(
                "data.size does not match adjacency. Expected %s"%self.ptr.data_size())
        cdef MapX3d* m = new MapX3d(&data[0,0],data.shape[0])
        self.ptr.init_data(m[0])
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
        cdef np.ndarray[np.double_t,ndim=2] res = np.zeros([self.ndata,3],dtype=np.double)
        cdef MapX3d* m = new MapX3d(&res[0,0],self.ndata)
        self.ptr.get_representer(m[0])
        del m
        return res

cdef class QGC_F3:
    """Implements GC::QuadHierarichalContraction<float,3> from quad.hpp"""
    cdef QGC_VAR_F3* ptr
    cdef int ndata
    def __cinit__(self, float max_var, int level):
        self.ptr = new QGC_VAR_F3(max_var, level)
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

