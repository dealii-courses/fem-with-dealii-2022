# Theory and Practice of Finite Element Methods

## Laboratory \#1

Goal of the laboratory: 
- set up a modern programming environment
- familiarse with git, VSCode, CMake, GoogleTest, and C++

Tasks:
- split the source file `pythagoras.cc` into two files, one containing
the `main` function, and the other one containing the tests
- move them to a `source` directory, and name them `main.cc` and `point_tests.cc`
- modify the `CMakeLists.txt` file to ensure the program still compiles and run
- Add three more tests that checks that the three-dimensional version of the `Point<dim>` class works as expected