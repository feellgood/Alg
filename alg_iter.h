#ifndef ALG_ITER_H
#define ALG_ITER_H

#include <iomanip>
#include <vector>
namespace alg {

  /**  The Iteration object calculates whether the solution has reached the
       desired accuracy, or whether the maximum number of iterations has
       been reached. This piece of code was copied and simplified from gmm++

       The method finished() checks the convergence.  The first()
       method is used to determine the first iteration of the loop.
  */

double norm(const std::vector<double> & X);


  class iteration {
  protected :
    double rhsn;       /* Right hand side norm.                            */
    size_t maxiter; /* Max. number of iterations.                       */
    int noise;         /* if noise > 0 iterations are printed.             */
    double resmax;     /* maximum residu.                                  */
    double resminreach, resadd;
    double diverged_res; /* Threshold beyond which the iterative           */
                       /* is considered to diverge.                        */
    size_t nit;     /* iteration number.                                */
    double res;        /* last computed residu.                            */
    bool written;

  public :

    iteration(double r = 1.0E-8, int noi = 0, size_t mit = (size_t)(-1),double div_res = 1E200)
      : rhsn(1.0), maxiter(mit), noise(noi), resmax(r), diverged_res(div_res)
    { nit = 0; res = 0.0; written = false; resminreach = 1E200; resadd = 0.0; }

    void  operator ++(int) {  nit++; written = false; resadd += res; }
    void  operator ++() { (*this)++; }

    bool first(void) { return nit == 0; }

    /* get/set the "noisyness" (verbosity) of the solvers */
    int get_noisy(void) const { return noise; }
    void set_noisy(int n) { noise = n; }
    void reduce_noisy(void) { if (noise > 0) noise--; }

    double get_resmax(void) const { return resmax; }
    void set_resmax(double r) { resmax = r; }

    double get_res() const { return res; }
    void enforce_converged(bool c = true)
    { if (c) res = double(0); else res = rhsn * resmax + double(1); }

    double get_diverged_residual(void) const { return diverged_res; }
    void set_diverged_residual(double r) { diverged_res = r; }

    size_t get_iteration(void) const { return nit; }
    void set_iteration(size_t i) { nit = i; }
    
    size_t get_maxiter(void) const { return maxiter; }
    void set_maxiter(size_t i) { maxiter = i; }

    double get_rhsnorm(void) const { return rhsn; }
    void set_rhsnorm(double r) { rhsn = r; }
    
    bool converged(void) {
      return !std::isnan(res) && res <= rhsn * resmax;
    }
    bool converged(double nr) { 
      res = std::fabs(nr);
      resminreach = std::min(resminreach, res);
      return converged();
    }
    bool converged(const std::vector<double> &v)
    { double norm2 =alg::norm(v); return converged(norm2); }
    
bool diverged(void)
	 { return std::isnan(res) || (nit>=maxiter) || (res>=rhsn*diverged_res && nit > 4); }

    bool diverged(double nr) {
      res = std::fabs(nr);
      resminreach = std::min(resminreach, res);
      return diverged();
    }

    bool finished(double nr) {
	if (noise > 0 && !written) {
        double a = (rhsn == 0) ? 1.0 : rhsn;
        converged(nr);
        std::cout << " iter " << std::setw(3) << nit << " residual " << std::setw(12) << std::fabs(nr) / a << std::endl;
        written = true;
      }
      return (converged(nr) || diverged(nr));
    }
    bool finished_vect(const std::vector<double> &v)
    { double norm2 =alg::norm(v); return finished(norm2); }

  };

}

#endif
