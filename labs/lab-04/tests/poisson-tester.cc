#include <gtest/gtest.h>

#include <fstream>

#include "poisson.h"

using namespace dealii;

// Test Fixture for Poisson problem
class PoissonTester : public ::testing::Test, public Poisson
{
public:
  PoissonTester() = default;
};


TEST_F(PoissonTester, MakeGrid)
{
  make_grid();
}