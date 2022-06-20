# 1. Numerical Solution of PDEs Using the Finite Element Method

This is shared course between the SISSA PhD track on Mathematical Analysis,
Modeling, and Applications (math.sissa.it) and the Master in High Performance
Computing (www.mhpc.it), on pratcical applications of the finite element method
using the deal.II (www.dealii.org) library (roughly 20h course).

The goal of the course is to provide to the students advanced numerical tools
for the solution of partial differential equations using the finite element
methods, and state-of-the-art programming and best-practices knowledge for the
implementation of their own codes.

The course enables a PhD or MHPC student working on numerical analysis of PDEs
to implement a state-of-the-art adaptive finite element code, that runs in
parallel, using modern C++ libraries. The implementation will be based on the
`deal.II` library (www.dealii.org).

What you will learn:
- How to use a modern C++ IDE, to build and debug your codes
- How to use a large FEM library to solve complex PDE problems
- How to properly document your code using Doxygen
- How to use a proper Git workflow to develop your applications
- How to leverage GitHub actions, google tests, and docker images to test and
  deploy your application
- How MPI parallelisation works in a real life FEM applications

## Main links

Course main page, with schedule and up to date information
- https://www.math.sissa.it/course/phd-course/fem-with-dealii-2022

Course slides, notes, materials, and codes:
- https://github.com/dealii-courses/fem-with-dealii-2022

Course recordings:
- https://www.youtube.com/playlist?list=PLcvf2raG3YsFISq9nymXe48YjD9dNwnHi

## Course program

A tentative detailed program is shown below 
(this will be updated during the course to reflect the actual course content)

1.  Tools and background
2.  Working with VisualStudio code, docker, and github
3.  Introduction to deal.II

4.  General structure of a deal.II code (step-1/step-2)

5.  Solving a Poisson problem using deal.II (step-3/step-4)

6.  Construction of Manufactured solutions
7.  Working on successively (uniformly) refined grids
8.  Studying the convergence rates of FEM codes (step-5)

9.  Hanging nodes constraints
10. Adaptive finite element methods in deal.II (step-6)
11. SOLVE-ESTIMATE-MARK-REFINE loop

12. Domain decomposition VS algebraic decomposition
13. Splitting workload: partitioning with space filling curves
14. Distributed memory parallelization (step-40)

15. Moving from serial Poisson to parallel distributed Poisson 
16. Vector valued problems
17. Block Preconditioning