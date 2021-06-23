#include "alg.h"
#include "alg_lu.h"

namespace alg
{

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
