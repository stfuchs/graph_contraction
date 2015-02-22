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
  template<typename T>
  struct VariancePolicy
  {
    struct VertexProps
    {
      VertexProps() {}
      VertexProps(T const& data) : sum(data), sum_sqr(data*data) {}

      inline T representer(size_t n) const { return sum/float(n); }

      T sum;
      T sum_sqr;
    };

    static inline float cost(VertexProps const& a, VertexProps const& b, size_t na, size_t nb)
    {
      float n_inv = 1./float(na+nb);
      return n_inv*( (a.sum_sqr+b.sum_sqr) - n_inv*(a.sum+b.sum)*(a.sum+b.sum) );
    }

    static inline void merge(VertexProps& a, VertexProps const& b)
    {
      a.sum += b.sum;
      a.sum_sqr += b.sum_sqr;
    }
  };
}
#endif
