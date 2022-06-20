#include <gtest/gtest.h>

#include <fstream>

#include "step-3.h"

using namespace dealii;

// Test Fixture for step-3
class Step3Tester : public ::testing::Test, public Step3
{
public:
  Step3Tester() = default;
};


TEST_F(Step3Tester, MakeGrid)
{
  make_grid();
}