#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    auto print_only_once = false;

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    if (!print_only_once || (print_only_once && comm_rank == 0)){
        printf("Hello world (rank %d of %d)\n", comm_rank, comm_size);
    }

    MPI_Finalize();
}
