/**
 * @file   graph_contraction.hpp
 * @author Steffen Fuchs <steffen@steffen-ubuntu>
 * @date   Sat Feb 21 17:07:50 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef GC_GRAPH_CONTRACTION_HPP
#define GC_GRAPH_CONTRACTION_HPP

#include <vector>
#include <list>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <Eigen/Core>
#include <gc/policies.hpp>

// Notations:
//   templates:
//   - Vd: vertex_descriptor type
//   - V: vertex type
//   - G: graph type
//   - T: data type

namespace GC
{
  struct Vertex { size_t vid; std::list<size_t> ids; };
  struct Edge { size_t eid; };

  typedef boost::adjacency_list<boost::listS,boost::listS,boost::undirectedS,Vertex,Edge> GraphT;
  typedef typename boost::graph_traits<GraphT>::vertex_descriptor Vd;
  typedef typename boost::graph_traits<GraphT>::edge_descriptor Ed;

  template<typename Scalar>
  struct EdgeProps
  {
    Ed edge;
    Scalar cost;
    bool outdated;
    bool invalid;
  };

  template<typename Scalar, int Dim, template<typename,int> class Policy>
  struct GraphContraction
  {
    typedef Policy<Scalar,Dim> PolicyT;
    typedef typename PolicyT::VertexProps VertexPropsT;
    typedef EdgeProps<Scalar> EdgePropsT;
    
    struct Event { size_t eid; };
    struct EventCompare
    {
      std::vector<EdgePropsT>* p;
      bool operator()(Event const& a, Event const& b) { return (*p)[a.eid].cost > (*p)[b.eid].cost; }
    };

    GraphT g;
    std::vector<EdgePropsT> eprops;
    std::vector<VertexPropsT> vprops;
    std::vector<Event> que;

    EventCompare comp;
    Scalar max_cost_;

    GraphContraction(Scalar max_cost) : max_cost_(max_cost){}

    void init_grid_adjacency(int hh,int ww)
    {
      std::vector<std::pair<int,int> > edges;
      for (int h=0; h<hh-1; ++h) {
        for (int w=0; w<ww-1; ++w) {
          edges.push_back( {h*ww+w, h*ww+w+1} );
          edges.push_back( {h*ww+w, (h+1)*ww+w} );
        }
        edges.push_back( {(h+1)*ww-1, (h+2)*ww-1} );
      }
      for (int w=0; w<ww-1; ++w) {
        edges.push_back( {(hh-1)*ww+w, (hh-1)*ww+w+1} );
      }
      init_adjacency(edges,ww*hh);
    }

    void init_adjacency(std::vector<std::pair<int,int> > const& edges, int ndata)
    {
      vprops.resize(ndata);
      eprops.resize(edges.size());
      que.resize(edges.size());
      comp.p = &eprops;
      g = GraphT(edges.begin(),edges.end(), ndata);

      size_t id = 0;
      auto es = boost::edges(g);
      for (auto eit=es.first; eit!=es.second; ++eit, ++id) {
        g[*eit] = { id };
        eprops[id] = { *eit, 0, true, false };
        que[id] = { id };
      }
    }

    template<typename DataT>
    void init_data(DataT const& data)
    {
      size_t id = 0;
      auto vs = boost::vertices(g);
      for (auto vit=vs.first; vit!=vs.second; ++vit, ++id) {
        g[*vit] = { id, {id} };
        vprops[id] = VertexPropsT(data.row(id));
      }
    }

    void fit()
    {
      int iteration = 0;
      size_t eid;
      Vd u, v;
      while (!que.empty())
      {
        eid = que.front().eid;
        pop_que();
        if (eprops[eid].invalid) continue;

        boost::tie(u,v) = split(eprops[eid].edge);        
        if (eprops[eid].outdated)
        {
          eprops[eid].outdated = false;
          eprops[eid].cost = PolicyT::cost(vprops[g[u].vid],vprops[g[v].vid],
                                           g[u].ids.size(),g[v].ids.size());
          push_que(eid);
        }
        else if (eprops[eid].cost <= max_cost_)
        {
          if (boost::out_degree(u,g) > boost::out_degree(v,g))
            contract(u,v);
          else
            contract(v,u);
          ++iteration;
        }
      }
    }

    void contract(Vd const& a, Vd const& b)
    {
      auto oe = boost::out_edges(b,g);
      for (auto oe_it=oe.first; oe_it!=oe.second; ++oe_it) {
        size_t eid = g[*oe_it].eid;
        eprops[eid].outdated = true;
        Vd b_to = boost::target(*oe_it,g);
        if (b_to == a)
          eprops[eid].invalid = true;
        else
          eprops[eid].edge = boost::add_edge(b_to, a, {eid}, g).first;
      }
      PolicyT::merge(vprops[g[a].vid], vprops[g[b].vid]);
      g[a].ids.splice(g[a].ids.end(), g[b].ids);
      boost::clear_vertex(b,g);
      boost::remove_vertex(b,g);
    }

    template<typename OutT>
    void get_labels(OutT& out) const
    {
      //out.resize(vprops.size(),1);
      int c = 0;
      auto vs = boost::vertices(g);
      for (auto vit=vs.first; vit!=vs.second; ++vit) {
        for (auto it=g[*vit].ids.begin(); it!=g[*vit].ids.end(); ++it) {
          out[*it] = c;
        }
        ++c;
      }
    }
    template<typename OutT>
    void get_representer(OutT& out) const
    {
      //out.resize(vprops.size(),Dim);
      auto vs = boost::vertices(g);
      for (auto vit=vs.first; vit!=vs.second; ++vit)
      {
        auto rep = vprops[g[*vit].vid].representer(g[*vit].ids.size());
        for (auto it=g[*vit].ids.begin(); it!=g[*vit].ids.end(); ++it) {
          out.row(*it) = rep;
        }
      }
    }
  
    inline void push_que(size_t eid)
    {
      que.push_back( {eid} );
      std::push_heap(que.begin(), que.end(), comp);
    }

    inline void pop_que()
    {
      std::pop_heap(que.begin(), que.end(), comp);
      que.pop_back();
    }

    inline std::pair<Vd,Vd> split(Ed e) const
    {
      return std::make_pair(boost::source(e,g), boost::target(e,g));
    }

    inline int data_size() const { return vprops.size(); }
  };
}
#endif
