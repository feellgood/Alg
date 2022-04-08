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
 monitor over the successive iterations if an algorithm is converging or not
 */
  class iteration {
  protected :
 
	/** Right hand side norm. */
	double rhsn;      
    
/** Max. number of iterations. */
	size_t maxiter; 
    
/** maximum residu. */
	double resmax;     

/** iteration number. */
    size_t nit;     
   
/** last computed residu. */
    double res;        
    
/** true : info was written */
    bool written;

  public :
	/** constructor */
inline iteration(double r = 1.0E-8, bool noise = false, size_t mit = (size_t)(-1)) : rhsn(1.0), maxiter(mit), resmax(r), nit(0), res(0), written(false), verbose(noise) { }

/** if true iterations are printed. */
const bool verbose;

/** increment of the number of iterations */    
inline void operator ++(int) { nit++; written = false; }

/** operator increment */
inline void operator ++() { (*this)++; }

/** true if iterations are starting */
inline bool first(void) { return nit == 0; }

/** getter for resmax */
inline double get_resmax(void) const { return resmax; }

/** setter for resmax */
inline void set_resmax(double r) { resmax = r; }

/** getter for residu res */
inline double get_res() const { return res; }

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
inline bool converged(double nr) 
	{ 
      res = std::fabs(nr);
      return converged();
	}

/** monitor the convergence through a vector */
inline bool converged(const std::vector<double> &v) { return converged( alg::norm(v) ); }


/** returns true if the algo has converged according the convergence criterias through a norm value nr */
inline bool finished(double nr) {
	if (verbose && !written) {
        double a = (rhsn == 0) ? 1.0 : rhsn;
        converged(nr);
        std::cout << " iter " << std::setw(3) << nit << " residual " << std::setw(12) << std::fabs(nr) / a << std::endl;
        written = true;
      }
    return converged(nr);
    }
    
/** returns true if the algo has converged according the convergence criterias through a vector */
inline bool finished_vect(const std::vector<double> &v) { return finished( alg::norm(v) ); }

  };
}

#endif
