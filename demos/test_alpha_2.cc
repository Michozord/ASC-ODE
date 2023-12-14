#define _USE_MATH_DEFINES 
#include <cmath>          //has to be the FIRST include, otherwise does not work!
#include <nonlinfunc.h>
#include <ode.h>


using namespace ASC_ode;

// the pendulum with a length constraint

// Lagrange = -mg*x(1) - mg*x(3) + x(4) * (x(0)^2 + x(1)^2 - 1) + x(5) * ((x(0) - x(2))^2 + (x(1)-x(3))^2 - 1)
// dLagrange
class dLagrange : public NonlinearFunction
{
  size_t DimX() const override { return 6; }
  size_t DimF() const override { return 6; }
  
  void Evaluate (ASC_bla::VectorView<double> x, ASC_bla::VectorView<double> f) const override
  {
    f(0) = 2*x(0)*x(4) + 2*x(5)*(x(0)-x(2));
    f(1) = - 1 + 2 * x(4) * x(1) + 2* x(5) * (x(1) - x(3));
    f(2) = - x(5) * 2 * (x(0) - x(2));
    f(3) = -1 - x(5)*2*(x(1)-x(3));
    f(4) = x(0)*x(0) + x(1)*x(1) - 1;
    f(5) = (x(0) - x(2))*(x(0) - x(2)) + (x(1)-x(3))*(x(1)-x(3)) - 1;    
  }
  void EvaluateDeriv (ASC_bla::VectorView<double> x, ASC_bla::MatrixView<double, ColMajor> df) const override
  {
    df(0,0) = 2*x(4) + 2*x(5);
    df(0,1) = 0;
    df(1,0) = 0;
    df(0,2) = -2*x(5);
    df(2,0) = -2*x(5);
    df(3,0) = 0;
    df(0,3) = 0;
    df(0,4) = 2*x(0);
    df(4,0) = 2*x(0);
    df(5,0) = 2*(x(0)-x(2));
    df(0,5) = 2*(x(0)-x(2));
    df(1,1) = 2*x(4) + 2*x(5);
    df(1,2) = 0;
    df(2,1) = 0;
    df(1,3) = -2*x(5);
    df(3,1) = -2*x(5);
    df(1,4) = 2*x(1);
    df(4,1) = 2*x(1);
    df(1,5) = 2*(x(1) - x(3));
    df(5,1) = 2*(x(1) - x(3));
    df(2,2) = x(5)*2;
    df(2,3) = 0;
    df(3,2) = 0;
    df(2,4) = 0;
    df(4,2) = 0;
    df(2,5) = -2*(x(0) - x(2));
    df(5,2) = -2*(x(0) - x(2));
    df(3,3) = 2*x(5);
    df(3,4) = 0;
    df(4,3) = 0;
    df(5,3) = -2*(x(1)-x(3));
    df(3,5) = -2*(x(1)-x(3));
    df(4,4) = 0;
    df(4,5) = 0;
    df(5,4) = 0;
    df(5,5) = 0;
  }
};


int main()
{
  double tend = 20*2*M_PI;
  double steps = 1000;
  Vector<double> x { 1, 0, 1, -1, 0, 0 };
  Vector<double> dx { 0, 0, 0, 0, 0, 0 };
  Vector<double> ddx { 0, 0, 0, 0, 0, 0 };
  auto rhs = std::make_shared<dLagrange>();
  auto mass = std::make_shared<Projector>(6, 0, 4);
  std::ofstream ost;
  ost.open ("../py_tests/output_alpha_2.txt");
  SolveODE_Alpha (tend, steps, 0.8, x, dx, ddx, rhs, mass, 
                   // [](double t, VectorView<double> x) { cout << "t = " << t << ", x = " << x(0) << " " << x(1) << " " << x(2) << endl; }
                   [&ost](double t, VectorView<double> x) { ost << t << " " << x(0) << " " << x(1) << " " << x(2) <<  " " << x(3) <<"\n"; }                   
                   );
  ost.close();
}
