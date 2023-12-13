#ifndef ODE_h
#define ODE_h

#include <functional>
#include <exception>
//#include <calcinverse.hpp>

#include "Newton.h"


namespace ASC_ode
{
  
  // implicit Euler method for dy/dt = rhs(y)
  void SolveODE_IE(double tend, int steps,
                   VectorView<double> y, std::shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    auto yold = std::make_shared<ConstantFunction>(y);
    auto ynew = std::make_shared<IdentityFunction>(y.Size());
    auto equ = ynew-yold - dt * rhs;

    double t = 0;
    for (int i = 0; i < steps; i++)
      {
        NewtonSolver (equ, y);
        yold->Set(y);
        t += dt;
        if (callback) callback(t, y);
      }
  }

  // explicit Euler method for dy/dt = rhs(y)
  void SolveODE_EE(double tend, int steps,
                   VectorView<double> y, std::shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;

    double t = 0;
    for (int i = 0; i < steps; i++)
      {
        Vector<double> f(rhs->DimF());
        rhs->Evaluate(y, f);
        y = Vector<double>(y + dt * f);
        t += dt;
        if (callback) callback(t, y);
      }
  }

 // <template ASC_bla::ORDERING ORD>
  void SolveODE_RK(double tend, int steps,
                   VectorView<double> y, std::shared_ptr<NonlinearFunction> rhs,
                   Matrix<double, ColMajor> A, Vector<double> b,
                   std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    auto yold = std::make_shared<ConstantFunction>(y);
    int s = b.Size();   // number of stages
    int n = y.Size();
    auto k = std::make_shared<IdentityFunction>(y.Size() * s);
    std::shared_ptr<NonlinearFunction>* funs = new std::shared_ptr<NonlinearFunction> [s];
    for (size_t i = 0; i<s; i++){
      auto tmp = std::make_shared<BlockMatVec>(A, k, i);
      funs[i] = Compose(rhs, yold + dt * tmp);
    }
    auto block_f = std::make_shared<BlockFunction>(s, funs);
    auto equ = k - block_f;
    double t = 0;
    for (size_t i = 0; i < steps; i++)
      {
        Vector<double> k_0 (n * s);
        for(size_t j=0; j < n * s; j++){
          rhs->Evaluate(y, k_0.Range(j * n, (j+1) * n));
        }
        std::cout<<"1\n";
        NewtonSolver (equ, k_0);
        std::cout<<"2\n";
        Vector<double> incr(n);
        incr = 0.;
        for(size_t l=0; l<s; l++){
          incr = incr + b(l) * k_0.Range(l*n, (l+1)*n); 
        }
        y = y + dt * incr;
        yold->Set(y);
        t += dt;
        if (callback) callback(t, y);
      }
  }

  

  
  
  
  // Newmark and generalized alpha:
  // https://miaodi.github.io/finite%20element%20method/newmark-generalized/
  
  // Newmark method for  mass*d^2x/dt^2 = rhs
  void SolveODE_Newmark(double tend, int steps,
                        VectorView<double> x, VectorView<double> dx,
                        std::shared_ptr<NonlinearFunction> rhs,   
                        std::shared_ptr<NonlinearFunction> mass,  
                        std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    double gamma = 0.5;
    double beta = 0.25;

    Vector<double> a(x.Size());
    Vector<double> v(x.Size());

    auto xold = std::make_shared<ConstantFunction>(x);
    auto vold = std::make_shared<ConstantFunction>(dx);
    auto aold = std::make_shared<ConstantFunction>(x);
    rhs->Evaluate (xold->Get(), aold->Get());
    
    auto anew = std::make_shared<IdentityFunction>(a.Size());
    auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
    auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

    auto equ = Compose(mass, anew) - Compose(rhs, xnew);
    double t = 0;
    for (int i = 0; i < steps; i++)            
      {
        NewtonSolver (equ, a);
        xnew -> Evaluate (a, x);
        vnew -> Evaluate (a, v);
        xold->Set(x);
        vold->Set(v);
        aold->Set(a);
        t += dt;
        if (callback) callback(t, x);
      }
    dx = v;
  }




  // Generalized alpha method for M d^2x/dt^2 = rhs
  void SolveODE_Alpha (double tend, int steps, double rhoinf,
                       VectorView<double> x, VectorView<double> dx, VectorView<double> ddx,
                       std::shared_ptr<NonlinearFunction> rhs,   
                       std::shared_ptr<NonlinearFunction> mass,  
                       std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    double alpham = (2*rhoinf-1)/(rhoinf+1);
    double alphaf = rhoinf/(rhoinf+1);
    double gamma = 0.5-alpham+alphaf;
    double beta = 0.25 * (1-alpham+alphaf)*(1-alpham+alphaf);

    Vector<double> a(x.Size());
    Vector<double> v(x.Size());

    auto xold = std::make_shared<ConstantFunction>(x);
    auto vold = std::make_shared<ConstantFunction>(dx);
    auto aold = std::make_shared<ConstantFunction>(ddx);
    // rhs->Evaluate (xold->Get(), aold->Get()); // solve with M ???
    
    auto anew = std::make_shared<IdentityFunction>(a.Size());
    auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
    auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

    // auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - Compose(rhs, (1-alphaf)*xnew+alphaf*xold);
    auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - (1-alphaf)*Compose(rhs,xnew) - alphaf*Compose(rhs, xold);

    double t = 0;
    a = ddx;

    for (int i = 0; i < steps; i++)
      {
        NewtonSolver (equ, a);
        xnew -> Evaluate (a, x);
        vnew -> Evaluate (a, v);

        xold->Set(x);
        vold->Set(v);
        aold->Set(a);
        t += dt;
        if (callback) callback(t, x);
      }
    dx = v;
    ddx = a;
  }

  

}


#endif
