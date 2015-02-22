/**
 * @file   print.hpp
 * @author Steffen Fuchs <steffen@steffen-ubuntu>
 * @date   Sun Feb 22 18:12:52 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef GC_PRINT_HPP
#define GC_PRINT_HPP

#include <iostream>

namespace GC
{
  template<typename G>
  void print_graph(G const& g)
  {
    auto vs = boost::vertices(g);
    for (auto vit=vs.first; vit!=vs.second; ++vit)
    {
      std::cout << "V: " << g[*vit].vid << " -> ";
      auto oe = boost::out_edges(*vit, g);
      for (auto eit=oe.first; eit!=oe.second; ++eit)
      {
        std::cout << g[boost::target(*eit,g)].vid << ", ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }

  template<typename GC>
  void print_que(GC const& gc)
  {
    std::cout << "Que: " << gc.que.size() << std::endl;
    for (auto it=gc.que.begin(); it!=gc.que.end(); ++it)
    {
      auto v = gc.split(gc.eprops[it->eid].edge);
      std::cout << gc.eprops[it->eid].outdated <<", " 
                << gc.eprops[it->eid].invalid <<", "
                << gc.eprops[it->eid].cost
                << " (" << gc.g[v.first].vid <<","<< gc.g[v.second].vid <<")\n";
    }
  }
}
#endif
