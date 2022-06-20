/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * based on deal.II step-1
 */


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <cmath>
#include <fstream>
#include <iostream>

using namespace dealii;


//! Generate a hypercube, and output it as an svg file.
void
first_grid(Triangulation<2> &triangulation)
{
  GridGenerator::hyper_cube(triangulation);

  std::cout << "Number of original vertices: " << triangulation.n_vertices()
            << std::endl;

  triangulation.refine_global(4);

  std::cout << "Number of vertices after 4 refinements: "
            << triangulation.n_vertices() << std::endl;
  {
    std::ofstream out("grid-1.svg");
    GridOut       grid_out;
    grid_out.write_svg(triangulation, out);
    std::cout << "Grid written to grid-1.svg" << std::endl;
  }
  {
    std::ofstream out("grid-1.vtk");
    GridOut       grid_out;
    grid_out.write_vtk(triangulation, out);
    std::cout << "Grid written to grid-1.vtk" << std::endl;
  }
}


//! Generate a locally refined hyper_shell, and output it as an svg file.
void
second_grid(Triangulation<2> &triangulation)
{
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);

  // triangulation.reset_all_manifolds();

  for (unsigned int step = 0; step < 5; ++step)
    {
      std::ofstream out("grid-2-" + std::to_string(step) + ".vtk");
      GridOut       grid_out;
      grid_out.write_vtk(triangulation, out);

      for (auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_center =
                center.distance(cell->vertex(v));

              if (std::fabs(distance_from_center - inner_radius) <=
                  1e-6 * inner_radius)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }

      triangulation.execute_coarsening_and_refinement();
    }
}

//! Create an L-shaped domain with one global refinement, and write it on
// `third_grid.vtk`.  Refine the L-shaped mesh adaptively around the re-entrant
// corner three times (after the global refinement you already did), but with a
// twist: refine all cells with the distance between the center of the cell and
// re-entrant corner is smaller than 1/3.
void
third_grid(Triangulation<2> &)
{
  // Insert code here
}

//! Returns a tuple with number of levels, number of cells, number of active
// cells. Test this with all of  your meshes.
std::tuple<unsigned int, unsigned int, unsigned int>
get_info(const Triangulation<2> &)
{
  // Insert code here
  return std::make_tuple(0, 0, 0);
}

int
main()
{
  {
    Triangulation<2> triangulation;
    first_grid(triangulation);
  }
  {
    Triangulation<2> triangulation;
    second_grid(triangulation);
  }
}
