#ifndef MASS_SPRING_H
#define MASS_SPRING_H



#include <../src/nonlinfunc.h>
#include <../src/ode.h>

using namespace ASC_ode;

#include <../ASC-bla/src/vector.h>
using namespace ASC_bla;



template <int D>
class Mass
{
public:
  double mass;
  Vector<double> pos;
  Vector<double> vel = {0.0,0.0,0.0};
  Vector<double> acc = {0.0,0.0,0.0};
};


template <int D>
class Fix
{
public:
  Vector<double> pos;
};


class Connector
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};

std::ostream & operator<< (std::ostream & ost, const Connector & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  std::array<Connector,2> connections;
};

template <int D>
class MassSpringSystem
{
  std::vector<Fix<D>> fixes;
  std::vector<Mass<D>> masses;
  std::vector<Spring> springs;
  Vector<double> gravity=0.0;
public:
  void SetGravity (Vector<double> _gravity) { gravity = _gravity; }
  Vector<double> Gravity() const { return gravity; }
  
  Connector AddFix (Fix<D> p)
  {
    fixes.push_back(p);
    return { Connector::FIX, fixes.size()-1 };
  }

  Connector AddMass (Mass<D> m)
  {
    masses.push_back (m);
    return { Connector::MASS, masses.size()-1 };
  }
  
  size_t AddSpring (Spring s) // double length, double stiffness, Connector c1, Connector c2)
  {
    springs.push_back (s); // Spring{length, stiffness, { c1, c2 } });
    return springs.size()-1;
  }



  
  auto & Fixes() { return fixes; } 
  auto & Masses() { return masses; } 
  auto & Springs() { return springs; }

  void GetState (VectorView<double> values, VectorView<double> dvalues, VectorView<double> ddvalues)
  {
    auto valmat = values.AsMatrix(Masses().size(), D);
    auto dvalmat = dvalues.AsMatrix(Masses().size(), D);
    auto ddvalmat = ddvalues.AsMatrix(Masses().size(), D);    

    for (size_t i = 0; i < Masses().size(); i++)
      {
        valmat.Row(i) = Masses()[i].pos;
        dvalmat.Row(i) = Masses()[i].vel;
        ddvalmat.Row(i) = Masses()[i].acc;
      }
  }
  
  void SetState (VectorView<double> values, VectorView<double> dvalues, VectorView<double> ddvalues)
  {
    auto valmat = values.AsMatrix(Masses().size(), D);
    auto dvalmat = dvalues.AsMatrix(Masses().size(), D);
    auto ddvalmat = ddvalues.AsMatrix(Masses().size(), D);

    for (size_t i = 0; i < Masses().size(); i++)
      {
        Masses()[i].pos = valmat.Row(i);
        Masses()[i].vel = dvalmat.Row(i);
        Masses()[i].acc = ddvalmat.Row(i);        
      }
  }
};

template <int D>
std::ostream & operator<< (std::ostream & ost, MassSpringSystem<D> & mss)
{
  ost << "fixes:" << std::endl;
  for (auto f : mss.Fixes())
    ost << f.pos <<std::endl;

  ost << "masses: " << std::endl;
  for (auto m : mss.Masses())
    ost << "m = " << m.mass << ", pos = " << m.pos <<std::endl;

  ost << "springs: " <<std::endl;
  for (auto sp : mss.Springs())
    ost << "length = " << sp.length << "stiffness = " << sp.stiffness
        << ", C1 = " << sp.connections[0] << ", C2 = " << sp.connections[1] <<std::endl;
  return ost;
}


template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  virtual size_t DimX() const { return D*mss.Masses().size(); }
  virtual size_t DimF() const { return D*mss.Masses().size(); }
  
  virtual void Evaluate (VectorView<double> x, VectorView<double> f) const
  {
    f = 0.0;
    
    auto xmat = x.AsMatrix(mss.Masses().size(), D);
    auto fmat = f.AsMatrix(mss.Masses().size(), D);
    
    for (size_t i = 0; i < mss.Masses().size(); i++)
      fmat.Row(i) = mss.Masses()[i].mass*mss.Gravity();
    
    for (auto spring : mss.Springs())
      {
        auto [c1,c2] = spring.connections;
        Vector<double> p1 (xmat.Width()), p2(xmat.Width());
        if (c1.type == Connector::FIX)
          p1 = mss.Fixes()[c1.nr].pos;
        else
          p1 = xmat.Row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.Fixes()[c2.nr].pos;
        else
          p2 = xmat.Row(c2.nr);

        double force = spring.stiffness * (Vector<double>(p1 + (-1)*p2).L2Norm()-spring.length);
        Vector<double> dir12 = 1.0/(Vector<double>(p1 + (-1)*p2).L2Norm()) * (p2+(-1)*p1);
        if (c1.type == Connector::MASS)
          fmat.Row(c1.nr) = fmat.Row(c1.nr) + force*dir12;
        if (c2.type == Connector::MASS)
          fmat.Row(c2.nr) = fmat.Row(c2.nr) + (-1) * force*dir12;
      }

    for (size_t i = 0; i < mss.Masses().size(); i++)
      fmat.Row(i) = (1/ mss.Masses()[i].mass)*fmat.Row(i) ;
  }
  
  virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const
  {
    // TODO: exact differentiation
    double eps = 1e-8;
    Vector<double> xl(DimX()), xr(DimX()), fl(DimF()), fr(DimF());
    for (size_t i = 0; i < DimX(); i++)
      {
        xl = x;
        xl(i) -= eps;
        xr = x;
        xr(i) += eps;
        Evaluate (xl, fl);
        Evaluate (xr, fr);
        df.Col(i) = 1/(2*eps) * (fr + (-1) * fl);
      }
  }
  
};

#endif
