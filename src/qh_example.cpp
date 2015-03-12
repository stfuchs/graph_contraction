
#include <iomanip>
#include <chrono>
#include <gc/quad.hpp>

typedef double Scalar;
typedef Eigen::Matrix<Scalar,3,Eigen::Dynamic> Mat;
static const int dim=3;

int main(int argc, char** argv)
{
  using namespace std::chrono;
  //system_clock::time_point start, end;
  //duration<double> elapsed;
  
  int h=384;
  int w=512;
  Mat mat = Mat::Random(dim,h*w);
  Eigen::Map<Mat> map(&mat(0),dim,h*w);
  Mat res = Mat::Zero(3,h*w);
  Eigen::Map<Mat> mres(&res(0),dim,h*w);
  std::cout << "Data:";
  std::cout << std::fixed << std::setw(11) << std::setprecision(4) << " \n";

  auto start = system_clock::now();
  GC::QuadHierarchicalContraction<Scalar,dim> gc(.5,.5,0);
  gc.init_data(h,w,map);
  gc.fit();
  gc.get_representer(mres);
  auto end = system_clock::now();
  
  auto elapsed = duration_cast<milliseconds>(end-start);
  std::cout << "duration: " << elapsed.count() <<" ms"<<std::endl;
  return 0;
}
