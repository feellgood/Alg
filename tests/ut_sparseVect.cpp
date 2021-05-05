#define BOOST_TEST_MODULE sparseVectTest

#include <boost/test/unit_test.hpp>

#include "alg.h" 

BOOST_AUTO_TEST_CASE(Stupid)
{
float x = 1.0;
BOOST_CHECK(x != 0.0f);
}

BOOST_AUTO_TEST_CASE(sparseVect)
{
alg::sparseVect v;

v.push_back(1,3.14);

BOOST_CHECK(v.isSorted() == false);
}


BOOST_AUTO_TEST_CASE(sparseVect_kill_zero)
{
alg::sparseVect v;

if(v.isEmpty()) {std::cout << "v is empty"  << std::endl;}
else {std::cout << "v is not empty"  << std::endl;}
BOOST_CHECK(v.isEmpty());

v.push_back(1,0);
v.kill_zero();
std::cout<<"after kill"<<std::endl;
if(v.isEmpty()) {std::cout << "v is empty"  << std::endl;}
else {std::cout << "v is not empty"  << std::endl;}
BOOST_CHECK(v.isEmpty());
}

