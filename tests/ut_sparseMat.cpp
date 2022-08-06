#define BOOST_TEST_MODULE sparseMatTest

#include <boost/test/unit_test.hpp>
#include <vector>
#include <iostream>

#include "alg_sparseMat.h"

BOOST_AUTO_TEST_SUITE(ut_mat)

BOOST_AUTO_TEST_CASE(Stupid)
{
float x = 1.0;
BOOST_CHECK(x != 0.0f);
}

BOOST_AUTO_TEST_CASE(w_sparseMat_build)
{
alg::w_sparseMat m(5);
m.push_back(1,1,3.14);
m.push_back(1,1,1000);
m.push_back(1,0,1);
m.push_back(3,1,3.14);
m.push_back(0,0,42);
m.push_back(2,2,456);
m.push_back(3,3,798);
m.push_back(4,4,159);
m.push_back(4,0,2.0);
BOOST_CHECK(m.getDim() == (size_t)5);
std::cout << "m=" << m << std::endl;

alg::r_sparseMat r_m(m);

std::cout << "r_m built done." << std::endl;
BOOST_CHECK(r_m.getDim() == m.getDim());
}

BOOST_AUTO_TEST_CASE(w_sparseMat1)
{
const size_t dim = 5;
alg::w_sparseMat m(dim);
m.push_back(1,1,3.14);
m.push_back(1,1,1000);
m.push_back(1,0,1);
m.push_back(3,1,3.14);
m.push_back(0,0,42);
m.push_back(2,2,456);
m.push_back(3,3,798);
m.push_back(4,4,159);
m.push_back(4,0,2.0);
std::cout << "m=" << m << std::endl;

alg::r_sparseMat r_m(m);

std::cout << "r_m is :" << std::endl;

r_m.print();

std::vector<double> x = {0.5,0,0,0,1};
std::vector<double> y(dim);

std::cout << "start mult" << std::endl;
r_m.mult(x,y); // y = r_m * x

for(unsigned int i=0;i<y.size();i++) {std::cout << i <<':'<<y[i] <<" ; ";}

BOOST_CHECK( y[0] == (double)21 );
BOOST_CHECK( y[1] == (double)0.5 );
BOOST_CHECK( y[2] == (double)0 );
BOOST_CHECK( y[3] == (double)0 );
BOOST_CHECK( y[4] == (double)160 );
}

BOOST_AUTO_TEST_SUITE_END()
