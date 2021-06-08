#include "alg.h"
#include <iostream>

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

double cg_ilu(alg::r_sparseMat& A, std::vector<double> & x, std::vector<double> & b, alg::iteration &iter) 
{
double rho, rho_1(0.0);
const size_t DIM = x.size();
std::vector<double> p(DIM),q(DIM),r(DIM),z(DIM);    

r_sparseMat LU = A;
ilu(LU);  // approximated LU decomposition   

iter.set_rhsnorm(alg::norm(b));
	
r.assign(b.begin(),b.end());// r = b;
std::vector<double> v_temp(DIM); 
alg::mult(A,x,v_temp);// v_temp = A x;
alg::sub(v_temp,r);// r -= v_temp; donc r = b - A x;

lu_solve(LU, r, z);   // z = LU \ r
rho = alg::dot(z,r);//rho = vect_sp(z, r);
p.assign(z.begin(),z.end());//copy(z, p);

while (!iter.finished_vect(r)) {
      if (!iter.first()) { 
 	     lu_solve(LU, r, z);   // z = LU \ r
	     rho = alg::dot(z,r);//rho = vect_sp(z, r);
	     alg::scaled(rho/rho_1,p); // p *= (rho/rho1)
		 alg::add(z,p);// p += z	; donc  p = z + (rho/rho_1)*p        
		 }
      alg::mult(A, p, q);
          
      if (alg::dot(q,p)==0) {std::cout << "CG solver with ILU abort" << std::endl; break;}	
	  double a=rho/alg::dot(q,p); //a = rho / vect_sp(q, p);	
	  alg::scaled_add(p, +a, x);  //add(scaled(p, +a), x);
	  alg::scaled_add(q, -a, r);  //add(scaled(q, -a), r);
      rho_1 = rho;
      ++iter;
      }   
return alg::norm(r)/alg::norm(b);
}

}//end namespace alg
