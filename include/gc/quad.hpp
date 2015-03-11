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
    typedef Eigen::Matrix<int,1,Eigen::Dynamic> LabelMat;
    typedef Eigen::Array<Scalar,Dim,Eigen::Dynamic> Arr;
    typedef Eigen::Array<Scalar,Eigen::Dynamic,Eigen::Dynamic> ArrY;

    Scalar max_cost_;
    int levels_;
    BaseContr g;

    QuadHierarchicalContraction(Scalar max_cost, int levels)
      : max_cost_(.5*max_cost), levels_(levels), g(max_cost) {}

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
      int ch = height;
      int cw = width;
      std::vector<Arr> sum, sum_sqr;
      std::vector<Eigen::Array<Scalar,1,Eigen::Dynamic> > cost;
      sum.push_back(data.array());
      sum_sqr.push_back(data.array()*data.array());
      for(int i=0;i<levels_;++i)
      {
        if (int(ch%2) != 0 || int(cw%2) != 0) {
          std::cout << "Max level: "<< i << std::endl; break;
        }
        sum.push_back(Arr());
        sum_sqr.push_back(Arr());
        quad_reduce(sum[i], ch, cw, sum.back());
        quad_reduce(sum_sqr[i], ch, cw, sum_sqr.back());
        ch /= 2;
        cw /= 2;
        Scalar n_inv = 1./pow(4.,Scalar(i+1));
        cost.push_back( (n_inv*(sum_sqr.back() - n_inv*sum.back()*sum.back())).colwise().sum() );
      }      

      Vd invalid = boost::add_vertex({0,{0} }, g.g);
      std::vector<Vd> vds(height*width,invalid);
      // down propagation of hierarchies
      g.vprops.reserve(width*height);
      for(int i=cost.size(); i>0; --i)
      {
        int si = 1<<i; // inner stride at level i
        int wi = width/si; // outer stride at level i
        for(int j=0; j<cost[i-1].cols(); ++j)
        {
          int li = (j/wi*width + j%wi)*si; // leave index
          if ( vds[li]==invalid && (cost[i-1][j] < max_cost_) )
          {
            g.vprops.push_back( { sum[i].col(j), sum_sqr[i].col(j) } );
            Vd vnew = boost::add_vertex(g.g);
            g.g[vnew].vid = g.vprops.size()-1;
            for (int r=li; r<li+si*width; r+=width) {
              for (int c=r; c<r+si; ++c) {
                vds[c] = vnew;
                g.g[vnew].ids.push_back(c);
              }
            }
          }
        }
      }
      // finalize lowest level
      for (size_t j=0; j<height*width; ++j)
      {
        if (vds[j] == invalid)
        {
          g.vprops.push_back( { sum[0].col(j), sum_sqr[0].col(j) } );
          vds[j] = boost::add_vertex( { g.vprops.size()-1, {j} }, g.g);
        }
      }
      // create horizontal edges
      std::cout << height*(width-1)+width*(height-1) << " ";
      g.eprops.reserve(height*(width-1)+width*(height-1));
      for (int h = 0; h<height; ++h)
      {
        for (int w = h*width; w<(h+1)*width-1; ++w)
        {
          if (vds[w] != vds[w+1])
          {
            auto res = boost::add_edge(vds[w], vds[w+1], g.g);
            if (res.second)
            {
              Ed enew = res.first;
              g.eprops.push_back( { enew, 0, false, false } );
              g.g[enew].eid = g.eprops.size()-1;
            }
          }
        }
      }
      // create vertical edges
      for (int w = 0; w<width; ++w)
      {
        for (int h = w; h<(height-1)*width+w; h+=width)
        {
          if (vds[h] != vds[h+width])
          {
            auto res = boost::add_edge(vds[h], vds[h+width], g.g);
            if (res.second)
            {
              Ed enew = res.first;
              g.eprops.push_back( { enew, 0, false, false } );
              g.g[enew].eid = g.eprops.size()-1;
            }
          }
        }
      }
      std::cout << g.eprops.size() << std::endl;
      boost::clear_vertex(invalid,g.g);
      boost::remove_vertex(invalid,g.g);
      g.make_que();
    }

    inline void fit() { g.fit(); };
    inline void get_labels(Eigen::Map<LabelMat>& out) const {
      g.get_labels(out);
    }
    inline void get_representer(Eigen::Map<Mat>& out) const {
      g.get_representer(out);
    }
  };
}

#endif
