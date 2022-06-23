#include "alg_bicg_ilu.h"

double alg::bicg_ilu(alg::sparseMat& A, std::vector<double> & x, const std::vector<double> & b, alg::iteration &iter) 
{
double rho_1(0.0), rho_2(0.0), alpha(0.0), beta(0.0), omega(0.0);
const size_t DIM = x.size();
if (b.size()!=DIM){std::cout << "rhs size mismatch" << std::endl; exit(1);}

std::vector<double> p(DIM), phat(DIM), shat(DIM), r(DIM), rt(DIM), s(DIM), t(DIM), v(DIM);    

sparseMat LU = A;
ilu(LU);  // approximated LU decomposition  

iter.set_rhsnorm(alg::norm(b));
	
r.assign(b.begin(),b.end());// r = b;
alg::mult(A, x, v);         // v = A x;
alg::sub(v, r);             // r -= v; donc r = b - A x;

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
      alg::mult(A, phat, v);   //  v = A phat;
     
	  alpha=rho_1/alg::dot(v, rt);    // alpha = rho_1 /(v'*rtilde);
      s.assign(r.begin(), r.end());   // s = r
	  alg::scaled_add(v, -alpha, s);  // s = s -alpha v; donc s = r -alpha v

      if (iter.finished_vect(s)){
	     alg::scaled_add(phat, alpha, x); // x = x + alpha phat
         break;
         }

      lu_solve(LU, s, shat);   // shat = LU \ s
      alg::mult(A, shat, t);               //  t = A shat;

      omega = alg::dot(t, s)/alg::dot(t,t); // omega = (t'* s) / (t'*t);
      alg::scaled_add(phat, alpha, x); // x = x + alpha phat;
      alg::scaled_add(shat, omega, x); // x = x + omega shat;

      alg::scaled(omega, t); // t *= omega
      r.assign(s.begin(), s.end());  // r = s
      alg::sub(t, r);                // r -= t	; donc  r = s - omega t

      rho_2 = rho_1;
      ++iter;
      }   
return iter.get_res()/iter.get_rhsnorm();
}
