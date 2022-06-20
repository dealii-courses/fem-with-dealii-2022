#  Lab 02 - Triangulation, DoFHandler, FiniteElement
## Numerical Solution of PDEs Using the Finite Element Method

**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

###  Setup

- You will need the following programs:
  - Visual Studio Code with the `Remote Development` plugin installed
  - Paraview (ver > 5.4 would be better)
-  Check that the `.devcontainer` directory contains a configuration that works
   foryou. When you open the root directory of the repository, it should prompt
   you toreopen in container. If not, try explictily with the command `F1 ->
   Reopen folder incontainer`
- If opening in container worked, you should be able to run a terminal (`CTRL +
  ~`) and see something like
```
To run a command as administrator (user "root"), use "sudo <command>".
See "man sudo_root" for details.
dealii@3f25621373bd:/workspaces/fem-with-dealii-2022$ 
```
- The folder `labs/lab-02` contains a `CMakeLists.txt` file. If you run 
  `code labs/lab-02` VSCode should open the `lab-02` directory, and provide
  commands to configure and build your application, at the bottom of the 
  terminal window


### Setup if you have ubuntu/debian and DO NOT want to use Docker

In this case you can copy the commands from the docker image that you can find here:

- <https://github.com/dealii/docker-files/blob/master/dependencies-focal/Dockerfile>

to install all the dependencies. After this, you can install the latest version 
of the library using the commands you find here:

- <https://github.com/dealii/dealii/blob/master/contrib/docker/Dockerfile>

### Lab-02: step-1

1.  Read the documentation of step-1 at
    -   <https://www.dealii.org/current/doxygen/deal.II/step_1.html>
2. Compile and run `step-1` inside the VSCode container and look at the output.
3. Understand what happens if you now call `first_grid(triangulation)` and
   `second_grid(triangulation)` one after the other. Fix the problem you get,
   **without** creating a second triangulation, or without going out of scope,
   and make sure your main function works correctly.
4. Create an image of an L-shape domain (add a function `third_grid` to step-1)
   with one global refinement, and output it as `third_grid.vtk` using the
   `GridOut::write_vtk()` function.
5. Refine the L-shaped mesh adaptively around the re-entrant corner three times
   (after the global refinement you already did), but with a twist: refine all
   cells with the distance between the center of the cell and re-entrant corner
   is smaller than 2/3.
6. Create a helper function that takes a reference to a Triangulation and
   returns a tuple with  the following information: number of levels, number of
   cells, number of active cells. Test this with all of your meshes.
7. Un-comment the `triangulation.reset_all_manifolds()` line in `second_grid()`.
   What happens now? Write your findings as a comment under the line in the
   file, and leave the command commented out after you experimented with it.
8. Bonus: Create a mesh that represents the surface of a torus and refine
   it 2 times globally. Output to vtk format and check the output. Note
   that your Triangulation needs to be of type ``Triangulation<2,3>``,
   which we will discuss later this week.
9. Bonus: Take a look at step-49 and read the included .msh file in your
   modified step-1 program.
10. Modify the file `gtest.cc` and enable all `DISABLED_` tests, by removing the
    `DISABLED_` prefix. Make sure you add tests for each of the points above.

### Lab-02: step-2

 1. Read the documentation of step-2 at
    -  <https://www.dealii.org/current/doxygen/deal.II/step_2.html>

 1. Run step-2. Look at the sparsity patterns in your browser.

 2. How does the pattern change if you increase the polynomial degree from 1 to
    2 orto 3?

 3. How does the pattern change if you use a globally refined (say 3 times)
    unitsquare?

 4. Are these patterns symmetric? Why/why not?

 5. How many entries per row in the sparsity pattern do you expect for a
    Q1 element (assuming four cells are around each vertex)? Check that
    this is true for the mesh in b) (look for `row_length(i)` and output
    them for each row). Can you construct a 2d mesh (without hanging
    nodes) that has a row with more entries?

 6. How many entries per row in the sparsity pattern are there for Q2 and
    Q3 elements, again assuming four cells around each vertex?

 7. Print all entries for row 42 for the original renumbered sparsity
    pattern.

 8. Bonus: Compute and output statistics like the number of
    unknowns, bandwidth of the sparsity pattern, average number of
    entries per row, and fill ratio.

 9. Design a test for each of the points above, and add it to the `gtest.cc`
    file (this will require some restructuring of `step-2`).  
