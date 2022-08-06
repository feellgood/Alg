#include "alg_bicg_ilu_dir.h"

double alg::bicg_ilu_dir(alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, const std::vector<double>& xd, const std::vector<size_t>& ld, alg::iteration &iter) 
{
double rho_1(0.0), rho_2(0.0), alpha(0.0), beta(0.0), omega(0.0);
const size_t DIM = x.size();
if (rhs.size()!=DIM){std::cout << "rhs size mismatch" << std::endl; exit(1);}

std::vector<double> p(DIM), phat(DIM), shat(DIM), r(DIM), rt(DIM), s(DIM), t(DIM), v(DIM), b(DIM);    
b.assign(rhs.begin(), rhs.end());// b = rhs;

A.mult(xd, v); // v = A * xd
alg::sub(v, b);      // b = b - A xd
zeroFill(ld,b);
iter.set_rhsnorm(alg::norm(b));

r_sparseMat LU = A;
ilu(LU);  // approximated LU decomposition  

iter.set_rhsnorm(alg::norm(b));
	
r.assign(b.begin(),b.end());// r = b;
A.mult(x, v);         // v = A x;
alg::sub(v, r);             // r -= v; donc r = b - A x;

zeroFill(ld,r);
rt.assign(r.begin(),r.end()); // copy(r, rt);
p.assign(r.begin(),r.end()); // copy(r, p );

while (!iter.finished_vect(r)) {

      rho_1 = alg::dot(rt,r);// rho_1 = vect_sp(rt, r);
      if (!iter.first()) { 
         beta = (rho_1 / rho_2) * (alpha / omega);
	     alg::scaled(omega, v); // v *= omega
		 alg::sub(v, p);        // p -= v	; donc  p = p - omega v  

	     alg::scaled(beta, p);  // p *= beta
		 alg::add(r, p);        // p += r	; donc  p = r + beta p       
		}

      lu_solve(LU, p, phat);   // phat = LU \ p
      zeroFill(ld,phat);
      A.mult(phat, v);   //  v = A phat;
      zeroFill(ld,v);
      
	alpha=rho_1/alg::dot(v, rt);    // alpha = rho_1 /(v'*rtilde);
      s.assign(r.begin(), r.end());   // s = r
	  alg::scaled_add(v, -alpha, s);  // s = s -alpha v; donc s = r -alpha v

      if (iter.finished_vect(s)){
	     alg::scaled_add(phat, alpha, x); // x = x + alpha phat
         break;
         }

      lu_solve(LU, s, shat);   // shat = LU \ s
      zeroFill(ld,shat);
      A.mult(shat, t);               //  t = A shat;
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
return iter.get_res()/iter.get_rhsnorm();
}

