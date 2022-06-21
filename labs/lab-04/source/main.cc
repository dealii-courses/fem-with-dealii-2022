#include "poisson.h"

int
main()
{
  deallog.depth_console(2);
  Poisson laplace_problem;
  laplace_problem.run();
  return 0;
}