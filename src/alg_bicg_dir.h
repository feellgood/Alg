#include "alg.h"

namespace alg
{

/** biconjugate gradient with dirichlet condtions and diagonal preconditionner */

double bicg_dir(alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, const std::vector<double>& xd, const std::vector<size_t>& ld, alg::iteration &iter) 
{
double rho_1(0.0), rho_2(0.0), alpha(0.0), beta(0.0), omega(0.0);
const size_t DIM = x.size();
if (rhs.size()!=DIM){std::cout << "rhs size mismatch" << std::endl; exit(1);}

std::vector<double> p(DIM), phat(DIM), shat(DIM), r(DIM), rt(DIM), s(DIM), t(DIM), v(DIM), diag_precond(DIM), b(DIM);    
b.assign(rhs.begin(),rhs.end());// b = rhs;

// le preconditionneur diagonal est une matrice diagonale contenant les inverses des coeffs de diag(A), ici on va stocker les coefficients dans un std::vector
for(unsigned int i=0;i<diag_precond.size();i++)
	{ diag_precond[i] = 1.0/A(i,i); }

alg::mult(A, xd, v);
alg::sub(v, b);      // b = b - A xd

std::for_each(ld.begin(),ld.end(),[&b,&diag_precond](const size_t _i){ b[_i] = 0.0; diag_precond[_i] = 0.0; });
iter.set_rhsnorm(alg::norm(b));
	
r.assign(b.begin(),b.end());// r = b;
alg::mult(A, x, v);         // v = A x;
alg::sub(v, r);             // r -= v; donc r = b - A x;

std::for_each(ld.begin(),ld.end(),[&r](const size_t _i){ r[_i] = 0.0; });
rt.assign(r.begin(),r.end()); // copy(r, rt);
p.assign(r.begin(),r.end()); // copy(r, p );

double tol=iter.get_resmax();
std::cout << "tol : " << tol << std::endl;

while (!iter.finished_vect(r)) {
//while (alg::norm(r) > tol * alg::norm(b)) {
      //std::cout << iter.get_iteration()  << "\t " << alg::norm(r) / alg::norm(b) << std::endl;

      rho_1 = alg::dot(rt,r);// rho_1 = vect_sp(rt, r);
      if (!iter.first()) { 
         beta = (rho_1 / rho_2) * (alpha / omega);
	     alg::scaled(omega, v); // v *= omega
		 alg::sub(v, p);        // p -= v	; donc  p = p - omega v  

	     alg::scaled(beta, p);  // p *= beta
		 alg::add(r, p);        // p += r	; donc  p = r + beta p       
		}

      alg::p_direct(diag_precond, p, phat); // phat = M p;
      alg::mult(A, phat, v);                //  v = A phat;
      std::for_each(ld.begin(),ld.end(),[&v](const size_t _i){ v[_i] = 0.0; }); // v[mask] = 0
     
	  alpha=rho_1/alg::dot(v, rt); // alpha = rho_1 /(v'*rtilde);
      s.assign(r.begin(), r.end());   // s = r
	  alg::scaled_add(v, -alpha, s);  // s = s -alpha v; donc s = r -alpha v

      if (iter.finished_vect(s)){
     // if (alg::norm(s) < tol * alg::norm(b)) {
	     alg::scaled_add(phat, alpha, x); // x = x + alpha phat
         break;
         }

      alg::p_direct(diag_precond, s, shat);// shat = M s;
      alg::mult(A, shat, t);               //  t = A shat;
      std::for_each(ld.begin(),ld.end(),[&t](const size_t _i){ t[_i] = 0.0; }); // t[mask] = 0

      omega = alg::dot(t, s)/alg::dot(t,t); // omega = (t'* s) / (t'*t);
      alg::scaled_add(phat, alpha, x); // x = x + alpha phat;
      alg::scaled_add(shat, omega, x); // x = x + omega shat;

      alg::scaled(omega, t); // t *= omega
      r.assign(s.begin(), s.end());  // r = s
      alg::sub(t, r);                // r -= t	; donc  r = s - omega t

      rho_2 = rho_1;
      ++iter;
      }   
alg::add(xd, x); // x = x + xd
return alg::norm(r)/alg::norm(b);
}

}//end namespace alg
