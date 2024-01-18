#include "mass_spring.h"

int main()
{
  MassSpringSystem<2> mss;
  mss.SetGravity( {0,-9.81} );
  auto fA = mss.AddFix( { { 0.0, 0.0 } } );
  auto mA = mss.AddMass( { 1, { 1.0, 0.0 } } );
  mss.AddSpring ( { 1, 10, { fA, mA } }  );

  auto mB = mss.AddMass( { 1, { 2.0, 0.0 } } );
  mss.AddSpring ( { 1, 20, { mA, mB } } );
  
  std::cout << "mss: " << std::endl << mss << std::endl;


  double tend = 10;
  double steps = 1000;
  
  Vector<double> x(2*mss.Masses().size());
  Vector<double> dx(2*mss.Masses().size());  
  Vector<double> ddx(2*mss.Masses().size());  

  auto mss_func = std::make_shared<MSS_Function<2>> (mss);
  auto mass = std::make_shared<IdentityFunction> (x.Size());      

  mss.GetState (x, dx, ddx);
  std::cout << "hello peeps"<<std::endl;
  SolveODE_Newmark(tend, steps, x, dx,  mss_func, mass,
                   [](double t, ASC_bla::VectorView<double> x) { std::cout << "t = " << t
                                                             << ", x = " << ASC_bla::Vector<double>(x) << std::endl; });
                                                             
  std::cout << "hello peeps"<<std::endl;
}
