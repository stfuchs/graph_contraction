#ifndef GC_QUAD_HPP
#define GC_QUAD_HPP

#include <gc/core.hpp>
#include <gc/policies.hpp>
#include <math.h>

namespace GC
{
  
  template<typename Scalar, int Dim>
  struct QuadHierarchicalContraction
  {
    typedef Variance<Scalar,Dim> Policy;
    typedef GraphContraction<Scalar,Dim,Variance> BaseContr;
    typedef typename BaseContr::VertexPropsT VertexProps;
    typedef typename BaseContr::EdgePropsT EdgeProps;
    typedef Eigen::Matrix<Scalar,Dim,Eigen::Dynamic> Mat;
    typedef Eigen::Array<Scalar,Dim,Eigen::Dynamic> Arr;
    typedef Eigen::Array<Scalar,Eigen::Dynamic,Eigen::Dynamic> ArrY;

    Scalar max_cost_;

    QuadHierarchicalContraction(Scalar max_cost) : max_cost_(max_cost) {}

    void quad_reduce(Arr& in, int ch, int cw, Arr& out)
    {
      int nw = cw/2;
      int nh = ch/2;

      typedef Eigen::Stride<Eigen::Dynamic,1> YStride;
      typedef Eigen::Map<Arr,0,Eigen::Stride<2*Dim,1> > XView;
      typedef Eigen::Map<ArrY,0,YStride> YView;
      
      Scalar* p = &in(0);
      XView x0 = XView(p,     Dim, nw*ch);
      XView x1 = XView(p+Dim, Dim, nw*ch);
      Arr sumx = x0+x1;
      p = &sumx(0);
      YView y0 = YView(p,        Dim*nw, nh, YStride(2*Dim*nw,1));
      YView y1 = YView(p+Dim*nw, Dim*nw, nh, YStride(2*Dim*nw,1));
      ArrY sumy = y0+y1;
      p = &sumy(0);
      out = Eigen::Map<Arr>(p,Dim,nw*nh);
    }

    void init_data(int height, int width, Mat const& data)
    {
      int hi=log2(height);
      int wi=log2(width);
      for(int i=1; i<log2(height); ++i){
        if( height % (2<<i) != 0) { hi=i-1; break; }
      }
      for(int i=1; i<log2(width); ++i){
        if( width % (2<<i) != 0) { wi=i-1; break; }
      }
      int level = std::min(hi,wi);
      std::cout << "Your hierarchy will have at most " << level << " levels\n";
      if (level!=0) {
      std::cout << "with a resolution of "
                << width/double(2<<level) << " x "
                << height/double(2<<level) << " at the top\n";
      }

      int ch = height;
      int cw = width;
      std::vector<Arr> sum, sum_sqr;
      std::vector<Eigen::Array<Scalar,1,Eigen::Dynamic> > cost;
      std::vector<Eigen::Array<bool,1,Eigen::Dynamic> > split;
      sum.push_back(data.array());
      sum_sqr.push_back(data.array()*data.array());
      for(int i=0;i<=level;++i)
      {
        sum.push_back(Arr());
        sum_sqr.push_back(Arr());
        quad_reduce(sum[i], ch, cw, sum.back());
        quad_reduce(sum_sqr[i], ch, cw, sum_sqr.back());
        ch /= 2;
        cw /= 2;
        Scalar n_inv = 1./pow(4.,Scalar(i+1));
        cost.push_back( (n_inv*(sum_sqr.back() - n_inv*sum.back()*sum.back())).colwise().sum() );
        split.push_back( cost.back()>max_cost_ );
        if (split.back().all())
        {
          std::cout << "break at level " << i << std::endl;
          break;
        }
      }
    }
  };
}

#endif
