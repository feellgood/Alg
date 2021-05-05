#include <random>

#include <boost/progress.hpp>
#include <boost/format.hpp>
#include <boost/tokenizer.hpp>

#include "gmm/gmm.h"

#include "alg.h"

using namespace gmm;
typedef gmm::row_matrix	<std::vector<double> >   write_matrix;
typedef gmm::row_matrix	<std::vector<double> >    read_matrix;

void cg_dir(const read_matrix& A, std::vector<double> & x, const std::vector<double> & b, const std::vector<size_t>& ld, iteration &iter) 
{
double rho, rho_1(0.0);
std::vector<double> p(x.size()),q(x.size()),r(x.size()),z(x.size());
std::vector<double> diag_precond(x.size());    

// le preconditionneur diagonal est une matrice diagonale contenant les inverses des coeffs de diag(A), ici on va stocker les coefficients dans un std::vector
for(unsigned int i=0;i<diag_precond.size();i++)
	{ diag_precond[i] = 1.0/A(i,i); }

iter.set_rhsnorm(alg::norm(b));
	
r.assign(b.begin(),b.end());// r = b;
std::vector<double> v_temp(x.size()); 
mult(A,x,v_temp);// v_temp = A x;
alg::dec(v_temp,r);// r -= v_temp; donc r = b - A x;

std::for_each(ld.begin(),ld.end(),[&r,&diag_precond](const size_t _i){ r[_i] = 0.0; diag_precond[_i] = 0.0; });

alg::p_direct(diag_precond,r,z);//mult(P, r, z);
rho = alg::p_scal(z,r);//rho = vect_sp(z, r);
p.assign(z.begin(),z.end());//copy(z, p);

while (!iter.finished_vect(r)) {
      if (!iter.first()) { 
 	        alg::p_direct(diag_precond,r,z);//mult(P, r, z);
	        rho = alg::p_scal(z,r);//rho = vect_sp(z, r);
	        alg::scaled(rho/rho_1,p); // p *= (rho/rho1)
		alg::inc(z,p);// p += z	; donc  p = z + (rho/rho_1)*p        
		}
	      mult(A, p, q);
          
	std::for_each(ld.begin(),ld.end(),[&q](size_t _i){q[_i] = 0.0; } );	      
	double a=rho/alg::p_scal(q,p); //a = rho / vect_sp(q, p);	
	alg::scaled_inc(p, +a, x); //add(scaled(p, +a), x);
	alg::scaled_inc(q, -a, r);//add(scaled(q, -a), r);
      rho_1 = rho;
      ++iter;
          }   
}


int main()
{
boost::timer time;
const int  VERBOSE = 0;
const int MAXITER = 5000;

const int NOD=101;

std::vector<size_t> ld;
std::vector<double> Vd(NOD, 0.);
std::vector<double> Lw(NOD);
std::vector<double> Xw(NOD,0.0);


std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_int_distribution<> distrib(0, 1000);

write_matrix Kw(NOD, NOD);
Kw(0,     0    )=+1.0;  Kw(0,     1    )=-1.0;
Kw(NOD-1, NOD-2)=-1.0;  Kw(NOD-1, NOD-1)=+1.0;

for (int n=1; n<NOD-1; ++n){
    Kw(n, n-1)=-1.0;  Kw(n, n  )=+2.0; Kw(n, n+1)=-1.0;
//Use `distrib` to transform the random unsigned int generated by gen 

/*
    double r= double(distrib(gen))/1000.;
    if (r>0.7){
       ld.push_back(n);
       Vd[n]=1.0;
       }
*/
}
ld.push_back(0);
Vd[0]=0.0;

ld.push_back(NOD-1);
Vd[NOD-1]=1.0;

// valeurs de dirichlet inserees dans le vecteur solution
std::for_each(ld.begin(),ld.end(),[&Xw,&Vd] (size_t _i){ Xw[_i] = Vd[_i]; } );

read_matrix  Kr(NOD, NOD);    gmm::copy(Kw, Kr);

std::vector<double> Lr(NOD,0.0);//read_vector  Lr(NOD);         
Lr.assign(Lw.begin(),Lw.end());//gmm::copy(Lw, Lr);

mult(Kr, Xw, Lw);

alg::scaled(Lw,-1.0,Lr);//add(scaled(Lw, -1.0), Lr);//Lr = -Lw

std::cout << boost::format("%5t solving %50T.") << std::flush;
time.restart();

gmm::iteration iter(1e-6);
iter.set_maxiter(MAXITER);
iter.set_noisy(VERBOSE);

gmm::clear(Xw);
cg_dir(Kr, Xw, Lr, ld, iter); // Conjugate gradient with dirichlet conditions and diagonal preconditionner
std::cout << "finished " << iter.get_iteration() << std::endl << "time elapsed : "<< time.elapsed() << endl;

for (int i=0; i<NOD; i+=10) { std::cout << i << "\t" << Xw[i]+Vd[i] << "\t" << Vd[i] << std::endl; }

return 0;
}


