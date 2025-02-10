#include <mpi.h>
#include <string>
#include <algorithm>
#include <cstring>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    auto print_only_once = false;

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    std::string message = "Hello World";
    char message_buffer[20];

    if(comm_rank == 0){
        std::copy_n(message.c_str(), strlen(message.c_str()) + 1, message_buffer);
        for(auto other_rank = 1; other_rank < comm_size; ++other_rank){
            MPI_Send(message_buffer, strlen(message_buffer) + 1, MPI_CHAR, other_rank, 27, MPI_COMM_WORLD);
        }
    }
    else{
        MPI_Status status;
        MPI_Recv(message_buffer, 20, MPI_CHAR, 0, 27, MPI_COMM_WORLD, &status);
    }

    if (!print_only_once || (print_only_once && comm_rank == 0)){
        const auto fmessage = std::string(message_buffer) + " from rank " + std::to_string(comm_rank) + "\n";
        printf("%s", fmessage.c_str());
    }

    MPI_Finalize();
}
