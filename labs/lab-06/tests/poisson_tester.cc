#include "poisson_tester.h"

#include <gtest/gtest.h>

#include <fstream>
#include <sstream>

using namespace dealii;

using PoissonTestTypes = ::testing::Types<std::integral_constant<int, 2>,
                                          std::integral_constant<int, 3>>;


using Poisson2DTester = PoissonTester<std::integral_constant<int, 2>>;

TYPED_TEST_CASE(PoissonTester, PoissonTestTypes);

TYPED_TEST(PoissonTester, MakeGrid)
{
  // Output dimension
  std::cout << "Working on dim " << TypeParam::value << std::endl;
  this->make_grid();
}