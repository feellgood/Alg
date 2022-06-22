#include "alg_bicg_dir.h"

double alg::bicg_dir(alg::sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, const std::vector<double>& xd, const std::vector<size_t>& ld, alg::iteration &iter) 
{
double rho_1(0.0), rho_2(0.0), alpha(0.0), beta(0.0), omega(0.0);
const size_t DIM = x.size();
if (rhs.size()!=DIM){std::cout << "rhs size mismatch" << std::endl; exit(1);}

std::vector<double> p(DIM), phat(DIM), shat(DIM), r(DIM), rt(DIM), s(DIM), t(DIM), v(DIM), diag_precond(DIM), b(DIM);    
b.assign(rhs.begin(),rhs.end());// b = rhs;

// le preconditionneur diagonal est une matrice diagonale contenant les inverses des coeffs de diag(A), ici on va stocker les coefficients dans un std::vector
A.buildDiagPrecond(diag_precond);

alg::mult(A, xd, v);
alg::sub(v, b);      // b = b - A xd

zeroFill(ld,b);
zeroFill(ld,diag_precond);
iter.set_rhsnorm(alg::norm(b));
	
r.assign(b.begin(),b.end());// r = b;
alg::mult(A, x, v);         // v = A x;
alg::sub(v, r);             // r -= v; donc r = b - A x;

zeroFill(ld,r);
rt.assign(r.begin(),r.end()); // copy(r, rt);
p.assign(r.begin(),r.end()); // copy(r, p );

double tol=iter.get_resmax();
std::cout << "tol : " << tol << std::endl;

while (!iter.finished_vect(r)) {
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
      zeroFill(ld,v);
      alpha=rho_1/alg::dot(v, rt); // alpha = rho_1 /(v'*rtilde);
      s.assign(r.begin(), r.end());   // s = r
	alg::scaled_add(v, -alpha, s);  // s = s -alpha v; donc s = r -alpha v

      if (iter.finished_vect(s)){
 	     alg::scaled_add(phat, alpha, x); // x = x + alpha phat
         break;
         }

      alg::p_direct(diag_precond, s, shat);// shat = M s;
      alg::mult(A, shat, t);               //  t = A shat;
	zeroFill(ld,t);
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
