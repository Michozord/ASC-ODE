#define _USE_MATH_DEFINES 
#include <cmath>          //has to be the FIRST include, otherwise does not work!
#include <nonlinfunc.h>
#include <ode.h>


using namespace ASC_ode;

// the pendulum with a length constraint

// Lagrange = -f*y + lam*(x*x+y*y-1)
// dLagrange
class dLagrange : public NonlinearFunction
{
  size_t DimX() const override { return 3; }
  size_t DimF() const override { return 3; }
  
  void Evaluate (ASC_bla::VectorView<double> x, ASC_bla::VectorView<double> f) const override
  {
    f(0) = 2*x(0)*x(2);
    f(1) = 2*x(1)*x(2) - 1;
    f(2) = x(0)*x(0)+x(1)*x(1)-1;
    
  }
  void EvaluateDeriv (ASC_bla::VectorView<double> x, ASC_bla::MatrixView<double, ColMajor> df) const override
  {
    df(0,0) = 2*x(2);
    df(0,1) = 0;
    df(0,2) = 2*x(0);

    df(1,0) = 0;
    df(1,1) = 2*x(2);
    df(1,2) = 2*x(1);

    df(2,0) = 2*x(0);
    df(2,1) = 2*x(1);
    df(2,2) = 0;
  }
};


int main()
{
  double tend = 20*2*M_PI;
  double steps = 1000;
  Vector<double> x { 1, 0, 0, };
  Vector<double> dx { 0, 0, 0 };
  Vector<double> ddx { 0, 0, 0 };
  auto rhs = std::make_shared<dLagrange>();
  auto mass = std::make_shared<Projector>(3, 0, 2);
  std::ofstream ost;
  ost.open ("../py_tests/output_alpha.txt");
  SolveODE_Alpha (tend, steps, 0.8, x, dx, ddx, rhs, mass, 
                   // [](double t, VectorView<double> x) { cout << "t = " << t << ", x = " << x(0) << " " << x(1) << " " << x(2) << endl; }
                   [&ost](double t, VectorView<double> x) { ost << t << " " << x(0) << " " << x(1) << " " << x(2) << "\n"; }                   
                   );
  ost.close();
}
