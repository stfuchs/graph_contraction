from libcpp.vector cimport vector
from libcpp.pair cimport pair
#cimport libeigenmap as lem
from libeigenmap cimport MapXi, MapXd, MapX3f, MapX3d
cimport numpy as np

cdef extern from "../include/gc/core.hpp":
    cdef cppclass GC_VAR_D1C "GC::GraphContraction<double,1,GC::Variance>":
        GC_VAR_D1C(double) except+
        void init_grid_adjacency(int,int)
        void init_adjacency(vector[pair[int,int]]&, int)
        void init_data(MapXd&)
        void fit()
        void get_labels(MapXi&)
        void get_representer(MapXd&)
    cdef cppclass GC_VAR_F3C "GC::GraphContraction<float,3,GC::Variance>":
        GC_VAR_F3C(float) except+
        void init_grid_adjacency(int,int)
        void init_adjacency(vector[pair[int,int]]&, int)
        void init_data(MapX3f&)
        void fit()
        void get_labels(MapXi&)
        void get_representer(MapX3f&)
    cdef cppclass GC_VAR_D3C "GC::GraphContraction<double,3,GC::Variance>":
        GC_VAR_D3C(double) except+
        void init_grid_adjacency(int,int)
        void init_adjacency(vector[pair[int,int]]&, int)
        void init_data(MapX3d&)
        void fit()
        void get_labels(MapXi&)
        void get_representer(MapX3d&)

cdef class GC_D1:
    cdef GC_VAR_D1C* ptr
    def __cinit__(self, double max_var):
        self.ptr = new GC_VAR_D1C(max_var)
    def __dealloc__(self):
        del self.ptr
    def init_grid_adjacency(self, int width, int height):
        self.ptr.init_grid_adjacency(width, height)
    def init_adjacency(self, edges, ndata):
        """ set graph adjacency manually
        edges: list of pairs of vertex indices:
               [(v_idx1,v_idx2),(idx...,idx...),...]
        ndata: number of data points
        """
        self.ptr.init_adjacency(edges, ndata)
    cpdef fit(self, object[np.double_t, ndim=1] data):
        #cdef MapXd* m = new MapXd(&data[0],data.size)
        cdef MapXd m
        m.reset(&data[0],data.size)
        #self.ptr.init_data(m)
        self.ptr.fit()
        self.ndata = data.size
    cpdef get_labels(self):
        cdef np.ndarray[np.int_t,ndim=1] res = np.zeros(self.ndata)
        cdef MapXi m
        #m.reset(<int*>res[0],self.ndata)
        #self.ptr.get_labels(m)
        return res
    cpdef get_representer(self):
        cdef np.ndarray[np.double_t,ndim=1] res = np.zeros(self.ndata)
        cdef MapXd m
        #m.reset(&res[0],self.ndata)
        #self.ptr.get_representer(m)
        return res


cdef class GC_D3:
    cdef GC_VAR_D3C* ptr
    def __cinit__(self, double max_var):
        self.ptr = new GC_VAR_D3C(max_var)
    def __dealloc__(self):
        del self.ptr
    def init_grid_adjacency(self, int width, int height):
        self.ptr.init_grid_adjacency(width, height)
    def init_adjacency(self, edges, ndata):
        """ set graph adjacency manually
        edges: list of pairs of vertex indices:
               [(v_idx1,v_idx2),(idx...,idx...),...]
        ndata: number of data points
        """
        self.ptr.init_adjacency(edges, ndata)
    cpdef fit(self, object[np.double_t, ndim=2] data):
        if data.shape[1] != 3:
            raise ValueError("This is the 3D implementation! Your dimension is %s"%data.shape[1])
        cdef MapX3d m
        #m.reset(&data[0,0],data.shape[0])
        #self.ptr.init_data(m)
        #self.ptr.fit()
        self.ndata = data.shape[0]
    cpdef get_labels(self):
        cdef np.ndarray[np.int_t,ndim=1] res = np.zeros(self.ndata)
        cdef MapXi m
        #m.reset(<int*>res[0],self.ndata)
        #self.ptr.get_labels(m)
        return res
    cpdef get_representer(self):
        cdef np.ndarray[np.double_t,ndim=2] res = np.zeros([self.ndata,3])
        cdef MapX3d m
        #m.reset(&res[0,0],self.ndata)
        #self.ptr.get_representer(m)
        return res
