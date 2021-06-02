#define BOOST_TEST_MODULE sparseMatTest

#include <boost/test/unit_test.hpp>
#include <vector>
#include <iostream>

#include "alg.h"

BOOST_AUTO_TEST_SUITE(ut_mat)

BOOST_AUTO_TEST_CASE(Stupid)
{
float x = 1.0;
BOOST_CHECK(x != 0.0f);
}


BOOST_AUTO_TEST_SUITE_END()
