#ifndef ALG_ITER_H
#define ALG_ITER_H

/** \file alg_iter.h
\brief iteration class from GMM, with some adaptations and simplifications.

The Iteration object calculates whether the solution has reached the
       desired accuracy, or whether the maximum number of iterations has
       been reached.

       The method finished() checks the convergence.  The first()
       method is used to determine the first iteration of the loop.
  */


#include <iomanip>
#include <vector>
namespace alg {

double norm(const std::vector<double> & X);

/**
 \class iteration
 monitor over the successive iterations if a value is converging or diverging
 */
  class iteration {
  protected :
 
	/** Right hand side norm. */
	double rhsn;      
    
/** Max. number of iterations. */
	size_t maxiter; 

/** if noise > 0 iterations are printed. */
    	int noise;
    
/** maximum residu. */
	double resmax;     
    
/** minimum value of the residu along iterations */	
	double resminreach;

/** Threshold beyond which the iterative is considered to diverge. */
    	double diverged_res; 

/** iteration number. */
    size_t nit;     
   
/** last computed residu. */
    double res;        
    
/** true : info was written */
    bool written;

  public :
	/** constructor */
inline iteration(double r = 1.0E-8, int noi = 0, size_t mit = (size_t)(-1),double div_res = 1E200)
      : rhsn(1.0), maxiter(mit), noise(noi), resmax(r), diverged_res(div_res)
    { nit = 0; res = 0.0; written = false; resminreach = 1E200; }

/** increment of the number of iterations */    
inline void operator ++(int) { nit++; written = false; }

/** operator increment */
inline void operator ++() { (*this)++; }

/** true if iterations are starting */
inline bool first(void) { return nit == 0; }

    /** get the "noisyness" (verbosity) of the solvers */
inline int get_noisy(void) const { return noise; }

/** set the "noisyness" (verbosity) of the solvers */
inline void set_noisy(int n) { noise = n; }
    
/** reduce noisiness */
inline void reduce_noisy(void) { if (noise > 0) noise--; }

/** getter for resmax */
inline double get_resmax(void) const { return resmax; }

/** setter for resmax */
inline void set_resmax(double r) { resmax = r; }

/** getter for residu res */
inline double get_res() const { return res; }

/** force status to be convergent */
inline void enforce_converged(bool c = true)
    { if (c) res = double(0); else res = rhsn * resmax + double(1); }

/** getter for diverged_res */
inline double get_diverged_residual(void) const
	{ return diverged_res; }

/** setter for diverged_res */
inline void set_diverged_residual(double r) { diverged_res = r; }

/** getter for number of iterations */
inline size_t get_iteration(void) const { return nit; }

/** setter for the number of iterations */
inline void set_iteration(size_t i) { nit = i; }
    
/** getter for the maximum number of iterations */
inline size_t get_maxiter(void) const { return maxiter; }

/** setter for the maximum number of iterations */
inline  void set_maxiter(size_t i) { maxiter = i; }

/** getter for the right hand side norm value */
inline double get_rhsnorm(void) const { return rhsn; }

/** setter for the right hand side norm */
    void set_rhsnorm(double r) { rhsn = r; }
    
/** return the monitored algo has converged or not according criteria fixed by right hand side norm rhsn */
inline bool converged(void)
	{ return !std::isnan(res) && res <= rhsn * resmax; }

/** monitor the convergence through a number */
inline bool converged(double nr) { 
      res = std::fabs(nr);
      resminreach = std::min(resminreach, res);
      return converged();
    }

/** monitor the convergence through a vector */
inline bool converged(const std::vector<double> &v)
    { double norm2 =alg::norm(v); return converged(norm2); }
    
/** return if the monitored algo has diverged */
inline bool diverged(void)
	 { return std::isnan(res) || (nit>=maxiter) || (res>=rhsn*diverged_res && nit > 4); }

/** monitor the divergence through a number */
inline bool diverged(double nr) {
      res = std::fabs(nr);
      resminreach = std::min(resminreach, res);
      return diverged();
    }

/** returns true if the algo has converged according the convergence criterias through a norm value nr */
inline bool finished(double nr) {
	if (noise > 0 && !written) {
        double a = (rhsn == 0) ? 1.0 : rhsn;
        converged(nr);
        std::cout << " iter " << std::setw(3) << nit << " residual " << std::setw(12) << std::fabs(nr) / a << std::endl;
        written = true;
      }
      return (converged(nr) || diverged(nr));
    }
    
/** returns true if the algo has converged according the convergence criterias through a vector */
inline bool finished_vect(const std::vector<double> &v)
    { double norm2 =alg::norm(v); return finished(norm2); }

  };
}

#endif
