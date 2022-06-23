// very simple MPI demo
// compile normally and run with
// mpirun -n X ./main

#include <mpi.h>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int  name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  std::cout << "Hi, I am process " << rank << " of " << size
            << " and I am running on " << processor_name << std::endl;



  if (rank == 0)
    {
      for (int i = 1; i < size; ++i)
        {
          double     value = 0.0;
          MPI_Status status;

          MPI_Recv(&value,         // void *buf
                   1,              // int count
                   MPI_DOUBLE,     // MPI_Datatype datatype,
                   MPI_ANY_SOURCE, // int source
                   0,              // int tag
                   MPI_COMM_WORLD,
                   &status); // MPI_Status *status

          std::cout << "I got value = " << value << " from "
                    << status.MPI_SOURCE << std::endl;
        }
    }
  else
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(3000));
      srand(rank);
      double my_value = (rand() % 1000) / 1000.0;

      MPI_Send(&my_value,  // void *buf
               1,          // int count
               MPI_DOUBLE, // MPI_Datatype datatype
               0,          // int dest
               0,          // int tag
               MPI_COMM_WORLD);
    }

  MPI_Finalize();
}
