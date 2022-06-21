#include "poisson_tester.h"

#include <gtest/gtest.h>

#include <fstream>
#include <sstream>

using namespace dealii;

using PoissonTestTypes = ::testing::Types<std::integral_constant<int, 1>,
                                          std::integral_constant<int, 2>,
                                          std::integral_constant<int, 3>>;


using Poisson2DTester = PoissonTester<std::integral_constant<int, 2>>;

TYPED_TEST_CASE(PoissonTester, PoissonTestTypes);

TYPED_TEST(PoissonTester, MakeGrid)
{
  // Output dimension
  std::cout << "Working on dim " << TypeParam::value << std::endl;
  this->make_grid();
}



// Test only two dimensional code
TEST_F(Poisson2DTester, TestLinear)
{
  std::stringstream str;

  str << "subsection Poisson<2>" << std::endl
      << "  set Dirichlet boundary condition expression = x" << std::endl
      << "  set Dirichlet boundary ids                  = 0" << std::endl
      << "  set Finite element degree                   = 1" << std::endl
      << "  set Forcing term expression                 = 0" << std::endl
      << "  set Grid generator arguments                = 0: 1: false"
      << std::endl
      << "  set Grid generator function                 = hyper_cube"
      << std::endl
      << "  set Neumann boundary condition expression   = 0" << std::endl
      << "  set Neumann boundary ids                    = " << std::endl
      << "  set Number of global refinements            = 4" << std::endl
      << "  set Number of refinement cycles             = 1" << std::endl
      << "  set Output filename                         = poisson" << std::endl
      << "  set Problem constants                       = pi:3.14" << std::endl
      << "end" << std::endl;

  parse_string(str.str());
  make_grid();
  setup_system();
  assemble_system();
  solve();

  auto tmp = solution;
  VectorTools::interpolate(dof_handler, dirichlet_boundary_condition, tmp);

  tmp -= solution;

  ASSERT_NEAR(tmp.l2_norm(), 0, 1e-10);
}

// Test only two dimensional code
TEST_F(Poisson2DTester, TestQuadratic)
{
  std::stringstream str;

  str << "subsection Poisson<2>" << std::endl
      << "  set Dirichlet boundary condition expression = x^2" << std::endl
      << "  set Dirichlet boundary ids                  = 0" << std::endl
      << "  set Finite element degree                   = 2" << std::endl
      << "  set Forcing term expression                 = -2" << std::endl
      << "  set Grid generator arguments                = 0: 1: false"
      << std::endl
      << "  set Grid generator function                 = hyper_cube"
      << std::endl
      << "  set Neumann boundary condition expression   = 0" << std::endl
      << "  set Neumann boundary ids                    = " << std::endl
      << "  set Number of global refinements            = 4" << std::endl
      << "  set Number of refinement cycles             = 1" << std::endl
      << "  set Output filename                         = quadratic"
      << std::endl
      << "  set Problem constants                       = pi:3.14" << std::endl
      << "end" << std::endl;

  parse_string(str.str());
  make_grid();
  setup_system();
  assemble_system();
  solve();

  auto tmp = solution;
  VectorTools::interpolate(dof_handler, dirichlet_boundary_condition, tmp);

  tmp -= solution;

  ASSERT_NEAR(tmp.l2_norm(), 0, 1e-10);
}


// Test only two dimensional code
TEST_F(Poisson2DTester, TestMixedBC1)
{
  std::stringstream str;

  str << "subsection Poisson<2>" << std::endl
      << "  set Dirichlet boundary condition expression = x^2" << std::endl
      << "  set Dirichlet boundary ids                  = 1,2,3" << std::endl
      << "  set Finite element degree                   = 2" << std::endl
      << "  set Forcing term expression                 = -2" << std::endl
      << "  set Grid generator arguments                = 0: 1: true"
      << std::endl
      << "  set Grid generator function                 = hyper_cube"
      << std::endl
      << "  set Neumann boundary condition expression   = 0" << std::endl
      << "  set Neumann boundary ids                    = 0" << std::endl
      << "  set Number of global refinements            = 4" << std::endl
      << "  set Number of refinement cycles             = 1" << std::endl
      << "  set Output filename                         = quadratic"
      << std::endl
      << "  set Problem constants                       = pi:3.14" << std::endl
      << "end" << std::endl;

  parse_string(str.str());
  make_grid();
  setup_system();
  assemble_system();
  solve();

  auto tmp = solution;
  VectorTools::interpolate(dof_handler, dirichlet_boundary_condition, tmp);

  tmp -= solution;

  ASSERT_NEAR(tmp.l2_norm(), 0, 1e-10);
}



// Test only two dimensional code
TEST_F(Poisson2DTester, TestLinearWithHangingNodes)
{
  std::stringstream str;

  str << "subsection Poisson<2>" << std::endl
      << "  set Dirichlet boundary condition expression = x" << std::endl
      << "  set Dirichlet boundary ids                  = 0" << std::endl
      << "  set Finite element degree                   = 1" << std::endl
      << "  set Forcing term expression                 = 0" << std::endl
      << "  set Grid generator arguments                = 0: 1: false"
      << std::endl
      << "  set Grid generator function                 = hyper_cube"
      << std::endl
      << "  set Neumann boundary condition expression   = 0" << std::endl
      << "  set Neumann boundary ids                    = " << std::endl
      << "  set Number of global refinements            = 4" << std::endl
      << "  set Number of refinement cycles             = 1" << std::endl
      << "  set Output filename                         = lin_with_handing"
      << std::endl
      << "  set Problem constants                       = pi:3.14" << std::endl
      << "end" << std::endl;

  parse_string(str.str());
  make_grid();

  for (unsigned int i = 0; i < 2; ++i)
    {
      for (const auto &cell : triangulation.active_cell_iterators())
        if (cell->center().square() <= .25)
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
    }


  setup_system();
  assemble_system();
  solve();
  output_results(0);

  auto tmp = solution;
  VectorTools::interpolate(dof_handler, dirichlet_boundary_condition, tmp);

  tmp -= solution;

  ASSERT_NEAR(tmp.l2_norm(), 0, 1e-10);
}