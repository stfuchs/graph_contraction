
#include <iomanip>
#include <gc/quad.hpp>

typedef double Scalar;
typedef Eigen::Matrix<Scalar,3,Eigen::Dynamic> Mat;
static const int dim=3;

int main(int argc, char** argv)
{
  int h=24;
  int w=32;
  Mat mat = Mat::Random(dim,h*w);
  Eigen::Map<Mat> map(&mat(0),dim,h*w);
  
  std::cout << "Data:";
  std::cout << std::fixed << std::setw(11) << std::setprecision(4) << " \n";

  GC::QuadHierarchicalContraction<Scalar,dim> gc(5.);
  gc.init_data(h,w,map);
  return 0;
}
