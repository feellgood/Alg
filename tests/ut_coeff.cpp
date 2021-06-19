#define BOOST_TEST_MODULE coeffTest

#include <boost/test/unit_test.hpp>
#include <vector>
#include <iostream>

#include "alg.h" 

BOOST_AUTO_TEST_SUITE(ut_alg)

BOOST_AUTO_TEST_CASE(affectation)
{
std::cout<< "check implicit operator=" <<std::endl;
alg::v_coeff a(42,3.14);
alg::v_coeff b = a;

BOOST_CHECK( a.getVal() == (double) 42.0 );
BOOST_CHECK( b.getVal() == a.getVal() );
a.setVal(25.0);
BOOST_CHECK( a.getVal() == (double) 25.0 );
BOOST_CHECK( b.getVal() == (double) 42.0 );
}

BOOST_AUTO_TEST_SUITE_END()
