#  Lab 06 - AFEM and distibuted memory parallelisation
## The finite element method using deal.II

**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

By the end of this laboratory, you will have modified your Poisson code to run
in parallel using distributed memory parallelization, and you will have some
knowledge of Task based parallelization

## Lab-06

1. Instrument your `Poisson` solver (either from lab-05 or lab-06) with a
`TimerOutput` class, and extract timing information about each method of the
`Poisson` class (see
https://www.dealii.org/current/doxygen/deal.II/classTimerOutput.html) using
scoped timers. Make a couple of test runs, and keep aside both the parameter
file you used to run the code, and the timing results.

2. Repeat the basic MPI commands from the file `mpihello.cc` and understand
how things work

3. Read carefully the documentation of `step-40`, and adapt your poisson solver
from lab-05 to become a parallel distributed solver. lab-06 contains a
(possible) solution (don't cheat!).

4. Similar to what is shown in the lecture, visualize the view of the mesh from
each individual processor using ``GridOut::write_vtk`` and the "global" mesh.
Use 3 MPI tasks for this, write a test that checks the reproducibility of the
code, by veryfing that if the mpi rank is 3, then 

5. Create a simple mesh (hyper_cube refined twice globally), run with two MPI
tasks and print locally owned, locally active, and locally relevant ``IndexSet``
for each task. Make a test that verify that these numbers are consistent acros
runs when gtest is run with mpirun.

6. Switch to release mode, decide on a global refinement level that takes in the
order of 30-60 seconds to solve, and study assembly and solve time with
1,2,4,8,12,16 MPI tasks. Which is the fastest, do the timings make sense based
on how many cores your machine has?

7. Play with the test problem by switching to 3d and changing the geometry to
something interesting. Your choice!
