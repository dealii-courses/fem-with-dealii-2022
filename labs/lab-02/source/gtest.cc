#include <deal.II/base/point.h>

#include <deal.II/grid/tria.h>

#include <gtest/gtest.h>

#include <fstream>

using namespace dealii;

// Declare functions in step-1 and step-2
void
first_grid(Triangulation<2> &);
void
second_grid(Triangulation<2> &);
void
third_grid(Triangulation<2> &);
std::tuple<unsigned int, unsigned int, unsigned int>
get_info(const Triangulation<2> &);



TEST(Step1, Mark1)
{
  Triangulation<2> tria;
  first_grid(tria);
  ASSERT_TRUE(std::ifstream("grid-1.svg"));
}


// Remove the DISABLED_ prefix from the following tests when you have done the
// corresponding exercise
TEST(Step1, DISABLED_Mark2)
{
  Triangulation<2> tria;
  second_grid(tria);
  ASSERT_TRUE(std::ifstream("grid-2.svg"));
}


TEST(Step1, DISABLED_Mark3)
{
  Triangulation<2> tria;
  third_grid(tria);
  ASSERT_TRUE(std::ifstream("grid-3.vtk"));
}



TEST(Step1, DISABLED_Mark4)
{
  Triangulation<2> tria;
  first_grid(tria);
  auto [levels, cells, active_cells] = get_info(tria);
  EXPECT_EQ(levels, 5u);
  EXPECT_EQ(cells, 341u);
  EXPECT_EQ(active_cells, 256u);
}



TEST(Step1, DISABLED_Mark5)
{
  Triangulation<2> tria;
  second_grid(tria);
  const auto [levels, cells, active_cells] = get_info(tria);
  EXPECT_EQ(levels, 6u);
  EXPECT_EQ(cells, 1250u);
  EXPECT_EQ(active_cells, 940u);
}



TEST(Step1, DISABLED_Mark6)
{
  Triangulation<2> tria;
  third_grid(tria);
  auto [levels, cells, active_cells] = get_info(tria);
  EXPECT_EQ(levels, 5u);
  EXPECT_EQ(cells, 351u);
  EXPECT_EQ(active_cells, 264u);
}



int
main(int argc, char *argv[])
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
