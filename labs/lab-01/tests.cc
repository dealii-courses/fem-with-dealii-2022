#include <deal.II/base/point.h>

#include <gtest/gtest.h>

using namespace dealii;


TEST(Pythagoras, Norm)
{
  Point<2> x(3, 4);
  ASSERT_EQ(x.norm(), 5);
}


TEST(Pythagoras, Distance)
{
  Point<2> x(4, 5);
  Point<2> y(1, 1);
  ASSERT_EQ(x.distance(y), 5);
}


TEST(Pythagoras, ScalarProduct)
{
  Point<2> x(3, 4);
  ASSERT_EQ(x * x, 25);
}

TEST(Pythagoras, Norm3D)
{
  Point<3> x(1, 2, 2);
  ASSERT_EQ(x.norm(), 3);
}


TEST(Pythagoras, Distance3D)
{
  Point<3> x(5, 1, 1);
  Point<3> y(1, 3, 5);
  ASSERT_EQ(x.distance(y), 6);
}


TEST(Pythagoras, ScalarProduct3D)
{
  Point<3> x(3, 4, 5);
  ASSERT_EQ(x * x, 50);
}