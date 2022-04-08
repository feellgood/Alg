#define BOOST_TEST_MODULE cg_dirTest

#include <boost/test/unit_test.hpp>
#include <vector>
#include <iostream>
#include <chrono>

#include "alg.h"

BOOST_AUTO_TEST_SUITE(ut_cg_dir)

BOOST_AUTO_TEST_CASE(cg_dir,* boost::unit_test::tolerance(1e-12))
{
auto t1 = std::chrono::high_resolution_clock::now();

const int  VERBOSE = 0;
const int MAXITER = 5000;
const int NOD=1001;

std::vector<size_t> ld;
std::vector<double> Vd(NOD, 0.);
std::vector<double> Lw(NOD);
std::vector<double> Xw(NOD,0.0);

std::cout<< "cg_dir: should find a linear solution" <<std::endl;

alg::w_sparseMat Kw(NOD);
Kw.push_back(0,0,1.0);
Kw.push_back(0,1,-1.0);
Kw.push_back(NOD-1,NOD-2,-1.0);
Kw.push_back(NOD-1,NOD-1,1.0);

for (int n=1; n<NOD-1; ++n)
	{
	Kw.push_back(n,n-1,-1.0);
	Kw.push_back(n,n,2.0);
	Kw.push_back(n,n+1,-1.0);
	}
ld.push_back(0);
Vd[0]=0.0;

ld.push_back(NOD-1);
Vd[NOD-1]=1.0;

// valeurs de dirichlet inserees dans le vecteur solution
//std::for_each(ld.begin(),ld.end(),[&Xw,&Vd] (size_t _i){ Xw[_i] = Vd[_i]; } );

alg::sparseMat Kr(Kw);

std::vector<double> Lr(NOD,0.0);
Lr.assign(Lw.begin(),Lw.end());

alg::mult(Kr,Xw,Lw);

alg::scaled(Lw,-1.0,Lr);//Lr = -Lw

alg::iteration iter(1e-6,VERBOSE,MAXITER);

Xw.clear();
Xw.resize(NOD);
double res = alg::cg_dir(Kr,Xw,Lr,Vd,ld,iter); // Conjugate gradient with dirichlet conditions and diagonal preconditionner

auto t2 = std::chrono::high_resolution_clock::now();

std::chrono::duration<double,std::micro> micros = t2-t1;
std::cout << "residu= " << res << "\tfinished " << iter.get_iteration() << std::endl << "time elapsed : " << micros.count() << " microsecondes\n";

for (int i=0; i<NOD; i+=50)
	{ 
	double val = Xw[i];//Xw[i]+Vd[i];
	double val_ref = (double) (i/((double) (NOD-1)));	
	//std::cout << i << "\t" << val << "\t" << Vd[i] << std::endl; 	
	std::cout << i << " : val = " << val << "\t" << "val_ref = " << val_ref << std::endl; 
	BOOST_TEST( val == val_ref );
	}
}

BOOST_AUTO_TEST_SUITE_END()
