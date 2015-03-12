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

#include <iostream>
#include <vector>
#include <list>
#include <chrono>
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
  template<typename Scalar, int Dim, template<typename,int> class Policy>
  struct GraphContraction
  {
    typedef Policy<Scalar,Dim> PolicyT;
    typedef typename PolicyT::Data DataT;
    typedef EdgeProps<Scalar> EdgePropsT;
    
    struct Event { size_t eid; };
    struct EventCompare
    {
      std::vector<EdgePropsT>* p;
      bool operator()(Event const& a, Event const& b) { return (*p)[a.eid].cost > (*p)[b.eid].cost; }
    };

    GraphT g;
    std::vector<EdgePropsT> eprops;
    std::vector<Event> que;

    EventCompare comp;
    Scalar max_cost;
    DataT* data;
    int edge_its;


    GraphContraction(Scalar _max_cost, DataT* _data)
      : max_cost(_max_cost), data(_data), edge_its(0){}

    void fit()
    {
      using namespace std::chrono;
      auto start = system_clock::now();
      int iteration = 0;
      size_t eid;
      Vd u, v;
      while (!que.empty())
      {
        eid = pop_que();
        if (eprops[eid].invalid) continue;

        boost::tie(u,v) = split(eprops[eid].edge);
        if (eprops[eid].outdated)
        {
          eprops[eid].cost = PolicyT::cost(g[u], g[v], data);
          if(eprops[eid].cost <= max_cost) {
            eprops[eid].outdated = false;
            push_que(eid);
          }
        }
        else if (eprops[eid].cost <= max_cost)
        {
          if (boost::out_degree(u,g) > boost::out_degree(v,g))
            contract(u,v);
          else
            contract(v,u);
          //++iteration;
        }
      }
      auto end = system_clock::now();
      auto elapsed = duration_cast<milliseconds>(end-start);
      std::cout << "fit duration   : " << elapsed.count() <<" ms"<<std::endl;
      /*std::cout << "collapsed edges: " << iteration << std::endl;
      std::cout << "moved edges    : " << edge_its << std::endl;
      std::cout << "num vertices   : " << boost::num_vertices(g) << std::endl;
      std::cout << "num edges      : " << boost::num_edges(g) << std::endl;
      std::cout << std::endl;*/
    }

    void contract(Vd const& a, Vd const& b)
    {
      //PolicyT::merge(g[a],g[b],data);
      //g[b].vid = g[a].vid;
      //return;
      
      auto oe = boost::out_edges(b,g);
      for (auto oe_it=oe.first; oe_it!=oe.second; ++oe_it) {
        //edge_its++;
        size_t eid = g[*oe_it].eid;
        eprops[eid].outdated = true;
        Vd b_to = boost::target(*oe_it,g);
        if (b_to == a)
          eprops[eid].invalid = true;
        else
          eprops[eid].edge = boost::add_edge(b_to, a, {eid}, g).first;
      }
      PolicyT::merge(g[a], g[b], data);
      boost::clear_vertex(b,g);
      boost::remove_vertex(b,g);
    }

    template<typename OutT>
    void get_labels(OutT& out) const
    {
      int c = 0;
      auto vs = boost::vertices(g);
      for (auto vit=vs.first; vit!=vs.second; ++vit) {
        for (auto it=g[*vit].ids.begin(); it!=g[*vit].ids.end(); ++it) {
          out(*it) = c;
        }
        ++c;
      }
    }
    template<typename OutT>
    void get_representer(OutT& out) const
    {
      auto vs = boost::vertices(g);
      for (auto vit=vs.first; vit!=vs.second; ++vit)
      {
        auto repr = PolicyT::repr(g[*vit], data);
        //for (auto it=data->vids[g[*vit].vid].begin(); it!=data->vids[g[*vit].vid].end(); ++it) {
        for (auto it=g[*vit].ids.begin(); it!=g[*vit].ids.end(); ++it) {
          out.col(*it) = repr;
        }
      }
    }

    inline void make_que()
    {
      Vd u,v;
      for (size_t i=0; i<eprops.size(); ++i)
      {
        boost::tie(u,v) = split(eprops[i].edge);
        eprops[i].cost = PolicyT::cost(g[u], g[v], data);
        if (eprops[i].cost <= max_cost){
          que.push_back({i});
          eprops[i].outdated = false;
        }
        else {
          boost::remove_edge(eprops[i].edge,g);
        }
      }
      comp.p=&eprops;
      std::make_heap(que.begin(),que.end(),comp);
    }
  
    inline void push_que(size_t eid)
    {
      que.push_back( {eid} );
      std::push_heap(que.begin(), que.end(), comp);
    }

    inline size_t pop_que()
    {
      size_t res = que.front().eid;
      std::pop_heap(que.begin(), que.end(), comp);
      que.pop_back();
      return res;
    }

    inline std::pair<Vd,Vd> split(Ed e) const
    {
      return std::make_pair(boost::source(e,g), boost::target(e,g));
    }
  };
}
#endif
