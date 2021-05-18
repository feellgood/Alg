#ifndef ALG_DENSEMAT_H
#define ALG_DENSEMAT_H

/** \file alg_denseMat.h 
 \brief dense matrix
 */

#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

namespace alg
{

/**
\class denseMat
usual dense matrix, container for values is a std::vector. Recommended for small dense matrices.

*/

  class denseMat : public std::vector<double> {
  protected:
    size_t _nrows, _ncols;

  public:

    inline const double& operator ()(size_t l, size_t c) const {
      assert(l < _nrows && c < _ncols);
      return *(this->begin() + c*_nrows+l);
    }
    inline double& operator ()(size_t l, size_t c) {
      assert(l < _nrows && c < _ncols);
      return *(this->begin() + c*_nrows+l);
    }

    const std::vector<double> &as_vector(void) const { return *this; }

    size_t&       nrows(void)        { return _nrows; }
    const size_t& nrows(void) const  { return _nrows; }

    size_t&       ncols(void)        { return _ncols; }
    const size_t& ncols(void) const  { return _ncols; }

    void swap(denseMat &m)
    { std::vector<double>::swap(m); std::swap(_ncols, m._ncols); std::swap(_nrows, m._nrows); }

    denseMat(std::vector<double> &v, size_t l, size_t c)
      : std::vector<double>(), _nrows(l), _ncols(c){
      assert(l*c == v.size());
      assign(v.begin(),v.end());
      }
 

    denseMat(size_t l, size_t c)
      : std::vector<double>(c*l), _nrows(l), _ncols(c)  { }
    denseMat(size_t l, size_t c, double value)
      : std::vector<double>(c*l, value), _nrows(l), _ncols(c)  { }
    denseMat(void) { _nrows = _ncols = 0; }

    /** clear function */
    inline void clear()
    { std::for_each(this->begin(),this->end(), [](double& x) { x=0.0; }); }

    /** printing function */
    inline void print(std::ostream & flux) const 
    { 
      for (size_t l=0; l<_nrows; ++l){
          flux<<'{'; 
          for (size_t c=0; c<_ncols; ++c){ 
              flux << *(this->begin() + c*_nrows+l) << " "; 
              }
          flux<<"}\n"; 
          }
    }
  };

/** operator<< for denseMat */
inline std::ostream & operator<<(std::ostream & flux, denseMat const& m) {m.print(flux); return flux;}
}
#endif //ALG_DENSEMAT_H

