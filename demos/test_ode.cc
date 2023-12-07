#include <nonlinfunc.h>
#include <ode.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
  size_t DimX() const override { return 2; }
  size_t DimF() const override { return 2; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -x(0);
  }
  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -1;
  }
};


int main()
{
  double tend = 4*M_PI;
  int steps = 1000;
  ASC_bla::Vector<double> y { 1, 0 };
  auto rhs = std::make_shared<MassSpring>();
  std::ofstream ost;
  ost.open ("C:/ESC/ASC-ODE/ASC-ODE/py_tests/output.txt");
  SolveODE_IE(tend, steps, y, rhs,
              [&ost](double t, VectorView<double> y) { ost << t << "  " << y(0) << " " << y(1) << "\n"; });
  ost.close();
}
