#  Lab 05 - Boundary conditions and constraints
## The finite element method using deal.II

**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

## General Instructions

For each of the point below, extend the `Poisson` class with functions that
perform the indicated tasks, trying to minimize the amount of code you copy and
paste, possibly restructuring existing code by adding arguments to existing
functions, and generating wrappers similar to the `run` method (e.g.,
`run_exercise_3`).

Once you created a function that performs the given task, add it to the
`poisson-tester.cc` file, and make sure all the exercises are run through the
`gtest` executable, e.g., adding a test for each exercise, as in the following
snippet:

```C++
TEST_F(PoissonTester, Exercise3) {
   run_exercise_3();
}
```

By the end of this laboratory, you will have modified your Poisson code to
mimic step-6 of the deal.II library, enabling local refinement.

## Lab-05

### step-5

1. Create a function that parses parameters from a string, to be used in the
testing infrastructure

2. Create a few tests that actually solve a Poisson problem on very small grids
with very simple but non-trivial combinations of boundary conditions, on
different domains, and verify the correctness of your findings. Examples
include using globally linear (quadratic) exact functions, with linear
(quadratic) finite elements, and verify that the error you make is actually
zero (in this case, global interpolation gives the exact solution, therefore
the finite element should also provide the exact solution)

3. Modify your code to use `AffineConstraints` instead of
`VectorTools::apply_boundary_values` (see the documentation of `step-6`). Run
again all tests, and verify that you pass your own checks again with the new
code

4. Add members functions
   
   void estimate(Vector<float> &error_per_cell) const;
   void mark_and_refine();

that implement the ESTIMATE-MARK-REFINE steps of tutorial step-6 of the deal.II
library in your poisson method. 

5. Add the parameters 

   double coarsening_fraction;
   double refinement_fraction;

and parse them from file (using the ParameterAcceptor infrastructure). If
refinement_fraction is 1.0, then global refinement should be used, and the
`estimate` function should not be called.

6. Add a preconditioner to your solver. If you have trilinos installed and it
is configured inside `deal.II`, use its algebraic multigrid preconditioner (it
works also with `deal.II` matrices). Otherwise use one of the other available
preconditioners in the library. Verify that your solver is now faster.

7. Plot a graph of the error VS number of degrees of freedom for a case with
non-trivial exact solution using both global refinement and local refinement.
What can you say about the order of convergence?

8. Make an estimate of the "cost per accuracy" in the two cases, i.e., plot the
computational cost required to reach a given accuracy in the two cases.