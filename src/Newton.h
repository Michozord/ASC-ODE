#ifndef Newton_h
#define Newton_h

#include "nonlinfunc.h"

namespace ASC_ode
{

  using namespace ASC_bla;

  void NewtonSolver (std::shared_ptr<NonlinearFunction> func, VectorView<double> x,
                     double tol = 1e-10, int maxsteps = 50,
                     std::function<void(int,double,VectorView<double>)> callback = nullptr)
  {
    Vector<double> res (func->DimF());
    Matrix<double, ColMajor> fprime(func->DimF(), func->DimX());

    for (int i = 0; i < maxsteps; i++)
      {
        func->Evaluate(x, res);
        std::cout << "|res| = " << std::endl;
        func->EvaluateDeriv(x, fprime);
        std::cout<<"is mir egal"<<std::endl;
        fprime = fprime.invert();
        //VectorView<double> tmp (fprime.Width(), );
        x = Vector<double>(x) + (-1)* Vector<double>(fprime*res);
        double err= res.L2Norm();
        if (callback)
          callback(i, err, x);
        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }

}

#endif
