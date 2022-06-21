#ifndef poisson_tester_h
#define poisson_tester_h

#include <gtest/gtest.h>

#include <fstream>

#include "poisson.h"

using namespace dealii;

// Test Fixture for Poisson problem, using integralconstant
template <class Integral>
class PoissonTester : public ::testing::Test, public Poisson<Integral::value>
{
public:
  PoissonTester() = default;
};

#endif