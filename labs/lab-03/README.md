#  Lab 03 - Poisson Problem
## The finite element method using deal.II

**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

## General Instructions

For each of the point below, extend the `step-3` class with functions that 
perform the indicated tasks, trying to minimize the amount of code you copy
and paste, possibly restructuring existing code by adding arguments to existing
functions, and generating wrappers similar to the `run` method (e.g., 
`run_exercise_3`).

Once you created a function that performs the given task, add it to the 
`step-3-tester.cc` file, and make sure all the exercises are run through
the `gtest` executable, e.g., adding a test for each exercise, as in the 
following snippet: 

```C++
TEST_F(Step3Tester, Exercise3) {
   run_exercise_3();
}
```

By the end of this laboratory, you will have a code that solves a Poisson 
problem in two dimensions, on different domain types, with arbitrary right hand
side, arbitrary Dirichlet boundary conditions, variable finite element degree, 
and different levels of refinement.

The program will read parameter files for its execution. This is the minimal 
starting point for simple but extensible Finite Element applications, and 
touches already many important and fundamental characteristics of general and 
extensible Finite Element codes.

## Lab-03 

### step-3

1.  See documentation of step-3 at
    <https://www.dealii.org/current/doxygen/deal.II/step_3.html>

2.  Compile and run step-3. Examine the source and header files. Move all 
    `#include<...>` directives that you can from the header `step-3.h` to
    the source file `step-3.cc`.

3.  Open the vtk output and visualize the solution in paraview or visit. 
    Figure out how to warp the solution by the solution variable. Save a picture
    of your visualization on the `figures` subdirectory, named `step-3-0.png`, 
    and commit it to your repository.

4.  Follow the instructions in “Modify the type of boundary condition”
    in the description of the tutorial, plot your output, and save a figure in
    the `figures`  subdirectory, named `step-3-1.png`, and commit it to your
    repository. Make sure the test is run through `gtest`.

5.  Now also do “A slight variation of the last point” but use the value
    -0.5 for the boundary with indicator 1. Save your output as `step-3-2.png`, 
    and commit it to the `figures` directory in your repository. Make sure the
    test is run through `gtest`.

6.  Change the setup to have $f=0$. Save your output as `step-3-3.png`, 
    and commit it to the `figures` directory in your repository. Make sure the
    test is run through `gtest`.

7.  Switch to an L-shaped domain and experiment with a combination of
    Dirichlet and Neumann boundary conditions. By experimentation, identify
    the faces adjacent to the re-entrant corner and apply Dirichlet conditions
    only there. Save your output as `step-3-4.png`, 
    and commit it to the `figures` directory in your repository. Make sure the
    test is run through `gtest`.

8.  Bonus: Do “Convergence of the mean”. Can you see the order $h^2$?
    Increase the polynomial order (you need to increase all orders of
    the quadratures in the program!) and check the convergence of the
    mean now. Make sure the test is run through `gtest`.

### ParameterHandler

1.  Follow the documentation of the class `ParameterAcceptor` 
    (from "A typical usage of this class"). Derive your `Step3` class from 
    `ParameterAcceptor`, and initialize the `ParameterAcceptor` class using 
    the string `"Step3"`. Add a function that initializes the class given a 
    parameter file name, i.e., `initialize("parameters.prm")`, by calling 
    `ParameterAcceptor::initialize(...)`
    
2.  Add the following parameters to your `Step3` class, and make sure they 
    are used throughout your code

      - Finite element degree
      - Number of global refinements
      - Output filename

    Notice that the `FE_Q` class will need to be built *after* the 
    initialization of the Finite element degree variable. You can create
    it right before you use it, and then reference to it by 
    `DoFHandler::get_fe()` when you need it later.

3.  Add two `FunctionParser<2>` members to your `Step3` class, one for the 
    boundary conditions, and one for the forcing term. Add two `std::string` 
    members, representing their mathematical expressions, and a 
    `std::map<std::string, double>` representing the constants to use in your
    mathematical expression, and make sure they are added as parameters of your 
    class as
    
      - Forcing term expression
      - Boundary condition expression
      - Problem constants
    
    and initialize the two `FunctionParser<2>` objects with the given 
    expressions and the given constants. Make sure your `Step3` class uses 
    these two functions where the boundary conditions are computed and where
    the right hand side is assembled

4.  Following the ideas in `step-70`, add two additional parameters to the 
    `Step3` class to select at run time what grid to create, and with what 
    arguments. Use as parameters

      - Grid generator function
      - Grid generator arguments

5.  Test your application with all the files in the `parameters` subdirectory,
    (running the test through the `gtest` application), and create a 
    visualization for each of the given parameters, with the same
    name of the parameter file, i.e., the visualization of the test that runs
    `parameters/hyper_shell.prm` should be saved on an image named `figures/hyper_shell.png`. Commits all your generated figures.
    Notice that for older versions of google test, you may need to run each
    test individually (this is a known issue between google testsuite and the
    `ParameterAcceptor` class). To run an individual test, you can use the 
    `--gtest_filter` command line option, i.e., the following command

    ```
    ./gtest --gtest_filter=Step3Tester.HyperShell
    ```

    will run only the test `Step3Tester.HyperShell`.

    If you want to run all tests one after the other you can call the `ctest` 
    utility, which will execute the `gtest` command passing the gtest filter for
    each of the tests you provided, producing, for example, an output similar 
    to:

    ```
    dealii@66c3be678611:/workspace/sissa-mhpc-lab-03/build-container$ ctest
    Test project /workspace/sissa-mhpc-lab-03/build-container
        Start 1: Step3Tester.MakeGrid
    1/6 Test #1: Step3Tester.MakeGrid .............   Passed    0.68 sec
        Start 2: Step3Tester.HyperCube
    2/6 Test #2: Step3Tester.HyperCube ............   Passed    0.78 sec
        Start 3: Step3Tester.HyperShell
    3/6 Test #3: Step3Tester.HyperShell ...........   Passed    1.90 sec
        Start 4: Step3Tester.HyperBall
    4/6 Test #4: Step3Tester.HyperBall ............   Passed    1.28 sec
        Start 5: Step3Tester.SinRhs
    5/6 Test #5: Step3Tester.SinRhs ...............   Passed    0.79 sec
        Start 6: Step3Tester.SinBc
    6/6 Test #6: Step3Tester.SinBc ................   Passed    1.26 sec
    
    100% tests passed, 0 tests failed out of 6
    ```