#  Lab 04 - Template parameters and convergence rates
## Theory and Practice of Finite Elements

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

By the end of this laboratory, you will have a code that solves a Poisson
problem in arbitrary dimensions, with Lagrangian finite elements of arbitrary
degree, on different domain types, with different boundary conditions, and
different functions for the definition of the right hand side, the stiffness
coefficient, and the forcing term.

The problem will run on successively refined grids, and we will verify
Bramble-Hilbert lemma for Lagrangian finite element spaces of different order,
building manufactured solutions using python, and plotting error convergence
tables using latex, tikz, and pgfplots.

The program will build on top of your implementation of Step3, drawing from
`step-4`, `step-5`, and `step-7`.

## Lab-04 

### step-4

1.  See documentation of step-4 at
    <https://www.dealii.org/current/doxygen/deal.II/step_4.html>

2.  Compile and run step-4. Examine the source and header files.

3. Copy your implementation of `step-3` from `lab-03` to the files
`source/poisson.cc`, `include/poisson.h`, and `tests/poisson-tester.cc`, make
sure you rename correctly all your files and classes to `Poisson`.

4. Add the template parameter `<int dim>` to your `Poisson` class, following
`step-4` as an example, and make sure that your program runs correctly both in
2D and in 3D.

5. Add the parameters
   
    - `Number of refinement cycles`
    - `Exact solution expression`
   
   and the corresponding member variables (i.e., `n_cycles`,
   `exact_solution_expression`, and `exact_solution`) and run the Poisson
   problem again for each refinement cycle with one global refinement, making
   sure you output the result for each refinement cycle separately in `vtu`
   format, i.e., if `Output filename` is `poisson_2d`, and `Number of
   refinement cycles` is 3, you should output

     - `poisson_2d_0.vtu`
     - `poisson_2d_1.vtu`
     - `poisson_2d_2.vtu`

    where the solution in `poisson_2d_0.vtu` should have `Number of global
    refinements` refinements, `poisson_2d_1.vtu` should have `Number of global
    refinements` +1 refinements, and `poisson_2d_2.vtu` should have `Number of
    global refinements` +2 refinements.

6. Add a `ParsedConvergenceTable` object to your `Poisson` class (see
https://www.dealii.org/current/doxygen/deal.II/classParsedConvergenceTable.html)
and add its parameters in the subsection `Error table` of the parameter file,
i.e., in the `Poisson` constructor add the following lines of code:
```
this->prm.enter_subsection("Error table");
error_table.add_parameters(this->prm);
this->prm.leave_subsection();
```

7. Set the boundary conditions, forcing function, and exact expression to get
the manufactured solution `u(x,y)=sin(pi*x)*cos(pi*y)`. Add a method
`compute_error()` to the `Poisson` class, that calls the
`ParsedConvergenceTable::error_from_exact` method with the `exact_solution`
function you created above. Make sure you output both the L2 and H1 error in
text format to a file. Play with the `jupyter` notebook
`manufactured_solutions.ipynb` to construct non-trivial exact solutions.

8. Add a parameter `Stiffness coefficient expression` and the corresponding
members to the `Poisson` class, so that the problem you will be solving is 
$-div(coefficient(x)\nabla u) = f(x)$.

9. Construct a (non-singular!) manufactured solution where `coefficient` is a
discontinuous function. Notice that the manufactured solution may need a
discontinuous forcing term on the right hand side, but should not have other
types of singularities, that is, you need to make sure that
$coefficient(x)\nabla u$ is continuous, i.e., that $\nabla u$ has a jump
depending on the jump of the `coefficient`. Output the error tables, and
comment on the error rates you observe. Do things improve when increasing the
polynomial order? Why?

10. (optional) Use the latex file provided in the latex subdirectory to
generate professional convergence plots for your solutions.