#include "alg_cg_ilu.h"

double alg::cg_ilu(alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & b, alg::iteration &iter) 
{
double rho, rho_1(0.0);
const size_t DIM = x.size();
std::vector<double> p(DIM),q(DIM),r(DIM),z(DIM);    

r_sparseMat LU = A; // ouch ...
ilu(LU);  // approximated LU decomposition   

iter.set_rhsnorm(alg::norm(b));
	
r.assign(b.begin(),b.end());// r = b;
std::vector<double> v_temp(DIM); 
A.mult(x,v_temp);// v_temp = A x;
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
      A.mult(p, q); // q = A * p
          
      if (alg::dot(q,p)==0) {std::cout << "CG solver with ILU abort" << std::endl; break;}	
	  double a=rho/alg::dot(q,p); //a = rho / vect_sp(q, p);	
	  alg::scaled_add(p, +a, x);  //add(scaled(p, +a), x);
	  alg::scaled_add(q, -a, r);  //add(scaled(q, -a), r);
      rho_1 = rho;
      ++iter;
      }   
return alg::norm(r)/alg::norm(b);
}
