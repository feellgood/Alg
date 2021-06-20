#include "alg.h"
#include "alg_sparseMat.h"

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
A.collect();  //  coefficients are sorted in lexicographic order

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


double cg_ilu_dir(alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, const std::vector<double> & xd, const std::vector<size_t>& ld, alg::iteration &iter) 
{
double rho, rho_1(0.0);
const size_t DIM = x.size();
if (rhs.size()!=DIM){std::cout << "rhs size mismatch" << std::endl; exit(1);}

std::vector<double> p(DIM),q(DIM),r(DIM),z(DIM),b(DIM);
b.assign(rhs.begin(),rhs.end());// b = rhs;

r_sparseMat LU = A;
ilu(LU);  // approximated LU decomposition    

alg::mult(A, xd, z); 
alg::sub(z, b);      // b = b - A xd

std::for_each(ld.begin(),ld.end(),[&b](const size_t _i){ b[_i] = 0.0; });
iter.set_rhsnorm(alg::norm(b));
	
r.assign(b.begin(),b.end());// r = b;
std::vector<double> v_temp(x.size()); 
alg::mult(A,x,v_temp);// v_temp = A x;
alg::sub(v_temp,r);   // r -= v_temp; donc r = b - A x;

std::for_each(ld.begin(),ld.end(),[&r](const size_t _i){ r[_i] = 0.0; });

lu_solve(LU, r, z);   // z = LU \ r
std::for_each(ld.begin(),ld.end(),[&z](const size_t _i){ z[_i] = 0.0; });

rho = alg::dot(z,r);        //rho = vect_sp(z, r);
p.assign(z.begin(),z.end());//copy(z, p);

while (!iter.finished_vect(r)) {
      if (!iter.first()) { 
         lu_solve(LU, r, z);   // z = LU \ r
         std::for_each(ld.begin(),ld.end(),[&z](const size_t _i){ z[_i] = 0.0; });

	     rho = alg::dot(z,r);
	     alg::scaled(rho/rho_1,p); // p *= (rho/rho1)
		 alg::add(z,p);            // p += z donc  p = z + (rho/rho_1)*p        
		 }
      alg::mult(A, p, q);
          
	  std::for_each(ld.begin(),ld.end(),[&q](size_t _i){q[_i] = 0.0; } );
      if (alg::dot(q,p)==0) {std::cout << "CG solver with ILU abort" << std::endl; break;}	      
	  double a=rho/alg::dot(q,p); //a = rho / vect_sp(q, p);	
	  alg::scaled_add(p, +a, x);  //add(scaled(p, +a), x);
	  alg::scaled_add(q, -a, r);  //add(scaled(q, -a), r);
      rho_1 = rho;
      ++iter; 
      }   
alg::add(xd, x); //x += xd
return alg::norm(r)/alg::norm(b);
}

}//end namespace alg
