#define _USE_MATH_DEFINES 
#include <cmath>        //has to be the FIRST include, otherwise does not work!
#include <nonlinfunc.h>
#include <ode.h>
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
    f(0) = x(1); // x(0) - position, x(1) - velocity
    f(1) = -x(0);
  }
  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -1;
  }
};

class Circuit : public NonlinearFunction{

  double resistance = 1;
  double capacity = 1;

  size_t DimX() const override {return 2;}
  size_t DimF() const override {return 2;}

  void Evaluate(VectorView<double> x, VectorView<double> f) const override 
  {
    //x(0) - time //x(1) - y(t)
    f(0)=1;   
    f(1)=( std::cos(100* M_PI* x(0))-x(1) ) / (resistance*capacity);
  }

  void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override
  {
    df = 0.0;
    df(1,0) =(-100 * M_PI* std::sin(100 * M_PI* x(0))-x(1))/(resistance* capacity) ;
    df(1,1) = -1 / (resistance*capacity);
  }
};


int main()
{
  double tend = 4*M_PI;
  int steps = 1000;
  ASC_bla::Vector<double> y { 1, 0 };
  auto rhs = std::make_shared<MassSpring>();
  std::ofstream ost;
  ost.open ("C:/Users/stein/Documents/ODE7.1.24/ASC-ODE/py_tests/output_ie.txt");
  SolveODE_IE(tend, steps, y, rhs,
              [&ost](double t, VectorView<double> y) { ost << t << "  " << y(0) << " " << y(1) << "\n"; });
  ost.close();
  y(0)=1; y(1)=0;
  std::ofstream ost2;
  ost2.open ("C:/Users/stein/Documents/ODE7.1.24/ASC-ODE/py_tests/output_ee.txt");
  SolveODE_EE(tend, steps, y, rhs,
              [&ost2](double t, VectorView<double> y) { ost2 << t << "  " << y(0) << " " << y(1) << "\n"; });
  ost2.close();
  y(0)=1; y(1)=0;
  std::ofstream ost3;
  ost3.open ("C:/Users/stein/Documents/ODE7.1.24/ASC-ODE/py_tests/output_cn.txt");
  SolveODE_CN(tend, steps, y, rhs,
              [&ost3](double t, VectorView<double> y) { ost3 << t << "  " << y(0) << " " << y(1) << "\n"; });
  ost3.close();

  y(0)=0; y(1)=0;
  tend=0.1;
  auto rhs_circuit = std::make_shared<Circuit>();
  std::ofstream ost4;
  ost4.open ("C:/Users/stein/Documents/ODE7.1.24/ASC-ODE/py_tests/output_circuit.txt");
  SolveODE_IE(tend, steps, y, rhs_circuit,
              [&ost4](double t, VectorView<double> y) { ost4 << t << "  " << y(1) << "\n"; });
  ost4.close();
}
