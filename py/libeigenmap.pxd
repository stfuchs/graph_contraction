cdef extern from "../src/eigenmap.cpp":
    cdef cppclass Map2i "EigenMap<int,2,1>":
        Map2i() except+
        Map2i(int*) except+
        void reset(int*)
    cdef cppclass Map2f "EigenMap<float,2,1>":
        Map2f() except+
        Map2f(float*) except+
        void reset(float*)
    cdef cppclass Map2d "EigenMap<double,2,1>":
        Map2d() except+
        Map2d(double*) except+
        void reset(double*)

    cdef cppclass Map3i "EigenMap<int,3,1>":
        Map3i() except+
        Map3i(int*) except+
        void reset(int*)
    cdef cppclass Map3f "EigenMap<float,3,1>":
        Map3f() except+
        Map3f(float*) except+
        void reset(float*)
    cdef cppclass Map3d "EigenMap<double,3,1>":
        Map3d() except+
        Map3d(double*) except+
        void reset(double*)
        
    cdef cppclass MapXi "EigenMap<int,Eigen::Dynamic,1>":
        MapXi() except+
        MapXi(int*, int) except+
        void reset(int*, int)
    cdef cppclass MapXf "EigenMap<float,Eigen::Dynamic,1>":
        MapXf() except+
        MapXf(float*, int) except+
        void reset(float*, int)
    cdef cppclass MapXd "EigenMap<double,Eigen::Dynamic,1>":
        MapXd() except+
        MapXd(double*, int) except+
        void reset(double*, int)
        
    cdef cppclass MapX2i "EigenMap<int,Eigen::Dynamic,2>":
        MapX2i() except+
        MapX2i(int*, int) except+
        void reset(int*, int)
    cdef cppclass MapX2f "EigenMap<float,Eigen::Dynamic,2>":
        MapX2f() except+
        MapX2f(float*, int) except+
        void reset(float*, int)
    cdef cppclass MapX2d "EigenMap<double,Eigen::Dynamic,2>":
        MapX2d() except+
        MapX2d(double*, int) except+
        void reset(double*, int)

    cdef cppclass MapX3i "EigenMap<int,Eigen::Dynamic,3>":
        MapX3i() except+
        MapX3i(int*, int) except+
        void reset(int*, int)
    cdef cppclass MapX3f "EigenMap<float,Eigen::Dynamic,3>":
        MapX3f() except+
        MapX3f(float*, int) except+
        void reset(float*, int)
    cdef cppclass MapX3d "EigenMap<double,Eigen::Dynamic,3>":
        MapX3d() except+
        MapX3d(double*, int) except+
        void reset(double*, int)

    cdef cppclass MapXXi "EigenMap<int,Eigen::Dynamic,Eigen::Dynamic>":
        MapXXi() except+
        MapXXi(int*, int, int) except+
        void reset(int*, int, int)
    cdef cppclass MapXXf "EigenMap<float,Eigen::Dynamic,Eigen::Dynamic>":
        MapXXf() except+
        MapXXf(float*, int, int) except+
        void reset(float*, int, int)
    cdef cppclass MapXXd "EigenMap<double,Eigen::Dynamic,Eigen::Dynamic>":
        MapXXd() except+
        MapXXd(double*, int, int) except+
        void reset(double*, int, int)
