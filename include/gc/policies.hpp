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
  struct Vertex { size_t vid; std::list<size_t> ids; };
  struct Edge { size_t eid; };

  typedef boost::adjacency_list<boost::listS,boost::setS,boost::undirectedS,Vertex,Edge> GraphT;
  typedef typename boost::graph_traits<GraphT>::vertex_descriptor Vd;
  typedef typename boost::graph_traits<GraphT>::edge_descriptor Ed;

  template<typename Scalar>
  struct EdgeProps
  {
    //EdgeProps() : edge(), cost(0), outdated(false), invalid(false) {}
    Ed edge;
    Scalar cost;
    bool outdated;
    bool invalid;
  };

  
  template<typename Scalar, int Dim>
  struct Variance
  {
    struct Data
    {
      std::vector<Scalar*> vsum;
      std::vector<Scalar*> vsumsqr;
      //std::vector<std::list<size_t> > vids;
    };
    
    typedef Eigen::Array<Scalar,Dim,1> Arr;
    typedef Eigen::Map<Arr> Map;

    static inline Scalar cost(Vertex const& a, Vertex const& b, Data* d)
    {
      Scalar n_inv = 1./Scalar(a.ids.size()+b.ids.size());
      Arr mean = Map(d->vsum[a.vid]) + Map(d->vsum[b.vid]);
      return ( n_inv*( Map(d->vsumsqr[a.vid]) + Map(d->vsumsqr[b.vid]) - n_inv*mean*mean ) ).sum();
    }

    static inline void merge(Vertex& a, Vertex& b, Data* d)
    {
      Map(d->vsum[a.vid]) += Map(d->vsum[b.vid]);
      Map(d->vsumsqr[a.vid]) += Map(d->vsumsqr[b.vid]);
      //d->vids[a.vid].splice(d->vids[a.vid].end(), d->vids[b.vid]);
      a.ids.splice(a.ids.end(), b.ids);
    }
    
    static inline Arr repr(Vertex const& a, Data* d)
    {
      //return Map(d->vsum[a.vid])/d->vids[a.vid].size();
      return Map(d->vsum[a.vid])/a.ids.size();
    }
  };
}
#endif
