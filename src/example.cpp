/**
 * @file   example.cpp
 * @author Steffen Fuchs <steffen@steffen-ubuntu>
 * @date   Sun Feb 22 17:00:41 2015
 * 
 * @brief  
 * 
 * 
 */

#include <iostream>
#include <iomanip>
#include <random>

#include "gc/core.hpp"
#include "gc/policies.hpp"



int main(int argc, char** argv)
{
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(0, 1);

  int ww = 10;
  int hh = 5;
  typedef Eigen::Matrix<double,Eigen::Dynamic,1> Mat;
  Mat data(ww*hh);
  std::vector<std::pair<size_t,size_t> > edges;

  std::cout << "Data:";
  std::cout << std::fixed << std::setw(11) << std::setprecision(4) << " \n";
  for (int h=0; h<hh; ++h) {
    for (int w=0; w<ww; ++w) {
      data[h*ww+w] = dist(mt);
      std::cout << data.coeff(h*ww+w) << "\t";
    }
    std::cout << "\n";
  }
  std::cout << "\n";


  Mat result(ww*hh,1);
  Eigen::Matrix<int,Eigen::Dynamic,1> ids(ww*hh,1);
  GC::GraphContraction<double, 1, GC::Variance> gc(.05);
  gc.init_grid_adjacency(hh,ww);
  gc.init_data(data);
  gc.fit();
  gc.get_representer(result);
  gc.get_labels(ids);

  for (int h=0; h<hh; ++h){
    for (int w=0; w<ww; ++w){
      std::cout << result.coeff(h*ww+w) << "\t";
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  for (int h=0; h<hh; ++h){
    for (int w=0; w<ww; ++w){
      std::cout << ids.coeff(h*ww+w) << "\t";
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  return 0;
}
