#include <iostream>
#include <Eigen/Core>


template<typename T, int R, int C>
struct EigenMap : Eigen::Map<Eigen::Matrix<T,R,C> >
{
  EigenMap()
    : Eigen::Map<Eigen::Matrix<T,R,C> >(NULL) { }
  EigenMap(T* data)
    : Eigen::Map<Eigen::Matrix<T,R,C> >(data) { }
  EigenMap(T* data, int cols)
    : Eigen::Map<Eigen::Matrix<T,R,Eigen::Dynamic> >(data, R, cols) { }
  EigenMap(T* data, int rows, int cols)
    : Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >(data, rows, cols) { }
  std::string rep() const
  {
    std::stringstream ss;
    ss << *this;
    return ss.str();
  }

  inline void reset(T* data) { new(this) EigenMap<T,R,C>(data); }
  inline void reset(T* data,int rows) { new(this) EigenMap<T,R,C>(data,rows); }
  inline void reset(T* data,int rows,int cols) { new(this) EigenMap<T,R,C>(data,rows,cols); }
};
