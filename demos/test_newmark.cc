#define _USE_MATH_DEFINES 
#include <cmath>        //has to be the FIRST include, otherwise does not work!
#include <nonlinfunc.h>
#include <ode.h>

using namespace ASC_ode;

class RHS : public NonlinearFunction
{
  size_t DimX() const override { return 1; }
  size_t DimF() const override { return 1; }
  
  void Evaluate (ASC_bla::VectorView<double> x, ASC_bla::VectorView<double> f) const override
  {
    f(0) = -x(0);
  }
  void EvaluateDeriv (ASC_bla::VectorView<double> x, ASC_bla::MatrixView<double, ColMajor> df) const override
  {
    df(0,0) = -1;
  }
};


int main()
{
  double tend = 2*M_PI;
  int steps = 100;
  Vector<double> x { 1, };
  Vector<double> dx { 0. };
  auto rhs = std::make_shared<RHS>();
  auto mass = std::make_shared<IdentityFunction>(1);
  SolveODE_Newmark(tend, steps, x, dx, rhs, mass,
                   [](double t, ASC_bla::VectorView<double> x) { std::cout << "t = " << t << ", x = " << x(0) << std::endl; }
                   );
}
