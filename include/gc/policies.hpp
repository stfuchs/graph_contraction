/**
 * @file   policies.hpp
 * @author Steffen Fuchs <steffen@steffen-ubuntu>
 * @date   Sun Feb 22 12:16:13 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef GC_POLICIES_HPP
#define GC_POLICIES_HPP

namespace GC
{
  template<typename Scalar, int Dim>
  struct Variance
  {
    typedef Eigen::Array<Scalar,Dim,1> Arr;
    typedef Eigen::Matrix<Scalar,Dim,1> Vec;
    typedef Eigen::Matrix<Scalar,Dim,Dim> Mat;

    struct VertexProps
    {
      VertexProps() {}
      VertexProps(Arr const& data)
        : sum(data), sum_sqr(data*data) {}
      VertexProps(Arr const& _sum, Arr const& _sum_sqr) :sum(_sum),sum_sqr(_sum_sqr) {} 

      inline Arr representer(size_t n) const { return (sum/Scalar(n)); }

      Arr sum;
      Arr sum_sqr;
    };

    static inline Scalar cost(VertexProps const& a, VertexProps const& b, size_t na, size_t nb)
    {
      Scalar n_inv = 1./Scalar(na+nb);
      Arr mean = a.sum+b.sum;
      return ( n_inv*((a.sum_sqr+b.sum_sqr) - n_inv*mean*mean) ).sum();
      //return cov.trace();
    }

    static inline void merge(VertexProps& a, VertexProps const& b)
    {
      a.sum += b.sum;
      a.sum_sqr += b.sum_sqr;
    }
  };

  template<typename Scalar>
  struct Variance<Scalar,1>
  {
    typedef Eigen::Matrix<Scalar,1,1> Vec;
    typedef Eigen::Matrix<Scalar,1,1> Mat;

    struct VertexProps
    {
      VertexProps() {}
      VertexProps(Vec const& data)
        : sum(data), sum_sqr(data*data) {}

      inline Vec representer(size_t n) const { return sum/Scalar(n); }

      Vec sum;
      Mat sum_sqr;
    };

    static inline Scalar cost(VertexProps const& a, VertexProps const& b, size_t na, size_t nb)
    {
      Scalar n_inv = 1./Scalar(na+nb);
      Vec mean = a.sum+b.sum;
      Mat cov = n_inv*((a.sum_sqr+b.sum_sqr) - n_inv*mean*mean);
      return cov.trace();
    }

    static inline void merge(VertexProps& a, VertexProps const& b)
    {
      a.sum += b.sum;
      a.sum_sqr += b.sum_sqr;
    }
  };
}
#endif
