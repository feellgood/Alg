#define BOOST_TEST_MODULE coreTest

#include <boost/test/unit_test.hpp>
#include <vector>
#include <iostream>

#include "alg.h"

BOOST_AUTO_TEST_SUITE(ut_core)

BOOST_AUTO_TEST_CASE(Stupid)
{
float x = 1.0;
BOOST_CHECK(x != 0.0f);
}

BOOST_AUTO_TEST_CASE(LinComb,* boost::unit_test::tolerance(1e-15))
{
const size_t DIM = 5;

alg::w_sparseMat A(DIM);
A.push_back(1,1,3.14);
A.push_back(1,1,1000);
A.push_back(1,0,1);
A.push_back(3,1,3.14);
A.push_back(0,0,42);
A.push_back(2,2,456);
A.push_back(3,3,798);
A.push_back(4,4,159);
A.push_back(4,0,2.0);
std::cout << "A=" << A << std::endl;

alg::sparseMat r_A(A);

std::vector<double> x {1,2,3,2,1};
std::vector<double> b {1,0,-3,1,2};

std::vector<double> v_temp(DIM),r(DIM);

r.assign(b.begin(),b.end());// r = b; 
alg::mult(r_A,x,v_temp);// v_temp = A x;
alg::sub(v_temp,r);// r -= v_temp; donc r = b - A x;

std::vector<double> r2(DIM);
alg::LinComb<false>(A,x,b,r2,std::minus<double>()); // r = b - A x
for(size_t i=0;i<DIM;i++) BOOST_TEST(r2[i] == r[i]);
}

BOOST_AUTO_TEST_SUITE_END()
