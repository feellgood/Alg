#include "alg.h"

namespace alg
{

inline std::ostream & operator<<(std::ostream & flux, std::vector<double> const& v) {
std::for_each(v.begin(),v.end(), [&flux](const double& x) { flux << x << " "; });
return flux;
}

void lu_solve(alg::r_sparseMat& LU, const std::vector<double> & b, std::vector<double> & x){
const size_t N=LU.getDim();
x.clear();
if (x.size() != N){ x.resize(N, 0.);}
 
// L y =b
std::vector<double> y(N, 0.0);
for (size_t i=0; i<N; i++) { 
    double scal=LU(i).dot(y);
    y[i]=b[i]-scal; // diag part of L matrix always equals 1
    }

// U x = Y
for (size_t i=N; i-- >0;) { // iterate from NOD-1 downto 0
    double scal=LU(i).dot(x);
    if (LU(i, i) == 0) {std::cout << "LU solver aborting" << std::endl; exit(1);}
    x[i]=(y[i]-scal)/LU(i, i);
    }
}

void ilu(alg::r_sparseMat& A){
size_t N=A.getDim();
for (size_t i=1; i<N; i++){ // row 0 unchanged  
    alg::sparseVect &Ai = A(i);
    double &Aii= Ai.getValRef(i);
    double Aii_bak=Aii;
    for (auto it_ik=Ai.begin(); it_ik != Ai.end() && it_ik->_i < i; ++it_ik){
        size_t k=it_ik->_i;
        double  Akk = A(k, k);
        double &Aik = it_ik->valRef(); // reference to A(i,k)
        if (Akk == 0) {std::cout << "Akk Nul" << std::endl; continue;}
        if (Aik == 0) {std::cout << "Aik Nul" << std::endl; continue;}
        Aik /= Akk;   

        auto it = std::find_if(Ai.begin(),Ai.end(),[&k](alg::v_coeff coeff){return (coeff._i > k); } );
        if (it==Ai.end()) {std::cout << "not found" << std::endl; break;}

        for (auto it_ij = it; it_ij != Ai.end(); ++it_ij){
            size_t j=it_ij->_i;
            double  &Aij = it_ij->valRef(); // reference to A(i,j)
            double  Akj = A(k,j);

            Aij -= Aik *Akj;    
            }

        }
    if (Aii==0) Aii=Aii_bak;
    
   }
}

void bicg_ilu(alg::r_sparseMat& A, std::vector<double> & x, std::vector<double> & b, alg::iteration &iter) 
{
double rho_1(0.0), rho_2(0.0), alpha(0.0), beta(0.0), omega(0.0);
const size_t DIM = x.size();
std::vector<double> p(DIM), phat(DIM), shat(DIM), r(DIM), rt(DIM), s(DIM), t(DIM), v(DIM);    

r_sparseMat LU = A;
ilu(LU);  // approximated LU decomposition  

iter.set_rhsnorm(alg::norm(b));
	
r.assign(b.begin(),b.end());// r = b;
alg::mult(A, x, v);         // v = A x;
alg::sub(v, r);             // r -= v; donc r = b - A x;

rt.assign(r.begin(),r.end()); // copy(r, rt);
p .assign(r.begin(),r.end()); // copy(r, p );

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
}

}//end namespace alg
