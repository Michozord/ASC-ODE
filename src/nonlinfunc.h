#ifndef NONLINFUNC_H
#define NONLINFUNC_H

#include <../ASC-bla/src/vector.h>
#include <../ASC-bla/src/matrix.h>


namespace ASC_ode
{
  using namespace ASC_bla;
  class NonlinearFunction
  {
  public:
    virtual ~NonlinearFunction() = default;
    virtual size_t DimX() const = 0;
    virtual size_t DimF() const = 0;
    virtual void Evaluate (VectorView<double> x, VectorView<double> f) const = 0;
    virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const = 0;
  };


  class IdentityFunction : public NonlinearFunction
  {
    size_t n;
  public:
    IdentityFunction (size_t _n) : n(_n) { } 
    size_t DimX() const override { return n; }
    size_t DimF() const override { return n; }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = x;
    }
    
    void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override
    {
      df = 0.0;
      df.Diag() = 1.0;
    }
  };



  class ConstantFunction : public NonlinearFunction
  {
    Vector<double> val;
    size_t dim_x;
  public:
    ConstantFunction (VectorView<double> _val) : val(_val), dim_x(val.Size()) { }
    ConstantFunction (VectorView<double> _val, size_t _dim_x) : val(_val), dim_x(_dim_x) { }
    void Set(VectorView<double> _val) { val = _val; }
    VectorView<double> Get() const { return val.View(); }
    size_t DimX() const override { return dim_x; }
    size_t DimF() const override { return val.Size(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = val;
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override
    {
      df = 0.0;
    }
  };

  
  
  class SumFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa, fb;
    double faca, facb;
  public:
    SumFunction (std::shared_ptr<NonlinearFunction> _fa,
                 std::shared_ptr<NonlinearFunction> _fb,
                 double _faca, double _facb)
      : fa(_fa), fb(_fb), faca(_faca), facb(_facb) { } 
    
    size_t DimX() const override { return fa->DimX(); }
    size_t DimF() const override { return fa->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      fa->Evaluate(x, f);
      f = faca * f;
      Vector<double> tmp(DimF());
      fb->Evaluate(x, tmp);
      f = f + facb*tmp;
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override
    {
      fa->EvaluateDeriv(x, df);
      Matrix<double, ColMajor> tmp(DimF(), DimX());
      tmp = faca * tmp;
      fb->EvaluateDeriv(x, tmp);
      df = df + facb*tmp;
    }
  };


  inline auto operator- (std::shared_ptr<NonlinearFunction> fa, std::shared_ptr<NonlinearFunction> fb)
  {
    return std::make_shared<SumFunction>(fa, fb, 1, -1);
  }

  inline auto operator+ (std::shared_ptr<NonlinearFunction> fa, std::shared_ptr<NonlinearFunction> fb)
  {
    return std::make_shared<SumFunction>(fa, fb, 1, 1);
  }

  
  class ScaleFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa;
    double fac;
  public:
    ScaleFunction (std::shared_ptr<NonlinearFunction> _fa,
                   double _fac)
      : fa(_fa), fac(_fac) { } 
    
    size_t DimX() const override { return fa->DimX(); }
    size_t DimF() const override { return fa->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      fa->Evaluate(x, f);
      f = fac*f;

    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override
    {
      fa->EvaluateDeriv(x, df);
      df = fac*df;
    }
  };

  inline auto operator* (double a, std::shared_ptr<NonlinearFunction> f)
  {
    return std::make_shared<ScaleFunction>(f, a);
  }




  // fa(fb)
  class ComposeFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa, fb;
  public:
    ComposeFunction (std::shared_ptr<NonlinearFunction> _fa,
                     std::shared_ptr<NonlinearFunction> _fb)
      : fa(_fa), fb(_fb) { } 
    
    size_t DimX() const override { return fb->DimX(); }
    size_t DimF() const override { return fa->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      Vector<double> tmp(fb->DimF());
      fb->Evaluate (x, tmp);
      fa->Evaluate (tmp, f);
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override
    {
      Vector<double> tmp(fb->DimF());
      fb->Evaluate (x, tmp);
      
      Matrix<double, ColMajor> jaca(fa->DimF(), fa->DimX());
      Matrix<double, ColMajor> jacb(fb->DimF(), fb->DimX());

      fb->EvaluateDeriv(x, jacb);
      fa->EvaluateDeriv(tmp, jaca);
      df = jaca*jacb;
    }
  };
  
  
  inline auto Compose (std::shared_ptr<NonlinearFunction> fa, std::shared_ptr<NonlinearFunction> fb)
  {
    return std::make_shared<ComposeFunction> (fa, fb);
  }
  
  class EmbedFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa;
    size_t firstx, dimx, firstf, dimf;
    size_t nextx, nextf;
  public:
    EmbedFunction (std::shared_ptr<NonlinearFunction> _fa,
                   size_t _firstx, size_t _dimx,
                   size_t _firstf, size_t _dimf)
      : fa(_fa),
        firstx(_firstx), dimx(_dimx), firstf(_firstf), dimf(_dimf),
        nextx(_firstx+_fa->DimX()), nextf(_firstf+_fa->DimF())
    { }
    
    size_t DimX() const override { return dimx; }
    size_t DimF() const override { return dimf; }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = 0.0;
      fa->Evaluate(x.Range(firstx, nextx), f.Range(firstf, nextf));
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override
    {
      df = 0;
      fa->EvaluateDeriv(x.Range(firstx, nextx),
                        df.Rows(firstf, nextf).Cols(firstx, nextx));      
    }
  };

  
  class Projector : public NonlinearFunction
  {
    size_t size, first, next;
  public:
    Projector (size_t _size, 
               size_t _first, size_t _next)
      : size(_size), first(_first), next(_next) { }
    
    size_t DimX() const override { return size; }
    size_t DimF() const override { return size; }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = 0.0;
      f.Range(first, next) = x.Range(first, next);
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double, ColMajor> df) const override
    {
      df = 0.0;
      df.Diag().Range(first, next) = 1; 
    }
  };

  class BlockFunction : public NonlinearFunction
  {
    size_t s;   // number of stages
    std::shared_ptr<NonlinearFunction>* funs;
  public:
    BlockFunction(size_t _s, std::shared_ptr<NonlinearFunction>* _funs)
      : s(_s), funs(_funs) { }

    size_t DimX() const override { return funs[0]->DimX(); }
    size_t DimF() const override { return s*funs[0]->DimF(); }

    void Evaluate(VectorView<double> x, VectorView<double> f) const override{
      for(size_t i=0; i<s; i++){
        std::shared_ptr<NonlinearFunction> fun = funs[i];
        size_t dim_f = fun->DimF();
        Vector<double> tmp(dim_f);
        fun->Evaluate(x, tmp);
        for(size_t j=0; j<dim_f; j++){
          f(i * dim_f + j) = tmp(j);
        }
      }
    }
    void EvaluateDeriv(VectorView<double> x, MatrixView<double, ColMajor> df) const override{
      size_t dim_f = funs[0]->DimF();
      Matrix<double, ColMajor> tmp (dim_f, s * dim_f);
      for(size_t j=0; j<s; j++){
        std::shared_ptr<NonlinearFunction> fun = funs[j];
        fun->EvaluateDeriv(x, tmp);     // tmp = Df_j Jacobi-Matrix of j-th function 
        df.Rows(j*dim_f, (j+1)*dim_f) = tmp;
      }
    }
  };

  class BlockMatVec : public NonlinearFunction
  {
    Matrix<double, ColMajor> A;
    std::shared_ptr<NonlinearFunction> vecfun;
    size_t j;
      
  public:
    BlockMatVec(Matrix<double, ColMajor> _A, std::shared_ptr<NonlinearFunction> _vecfun, size_t _j)
      : A(_A), vecfun(_vecfun), j(_j) { }
    
    size_t DimX() const override { return vecfun->DimX(); }
    size_t DimF() const override { return (vecfun->DimF())/A.Height(); }

    void Evaluate(VectorView<double> x, VectorView<double> f) const override{
      f = 0.;
      size_t s = A.Height();
      size_t n = (vecfun->DimF())/s;
      Vector<double> tmp (vecfun->DimF());
      vecfun->Evaluate(x, tmp);
      for(size_t l=0; l<s; l++){
        f = f + (A(j, l) * tmp.Range(n*l, n*(l+1)));
      }
    }
    void EvaluateDeriv(VectorView<double> x, MatrixView<double, ColMajor> df) const override {
      size_t s = A.Height();
      size_t n = (vecfun->DimF())/s;
      for(size_t l = 0; l < s; l++){
        Matrix<double, ColMajor> tmp (n, n);
        tmp.Diag() = A(j, l);
        df.Cols(n*l, n*(l+1)) = tmp;
      }
    }
  };

  
}

#endif
