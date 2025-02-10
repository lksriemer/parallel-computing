#include <algorithm>
#include <cstring>
#include <mpi.h>
#include <stdio.h>
#include <string>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  auto print_only_once_per = true;
  auto use_sendrecv = false;

  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  const std::string message =
      std::string("Hello World from thread ") + std::to_string(comm_rank);
  char send_buffer[100];
  char recv_buffer[100];
  std::fill(send_buffer, send_buffer + 100, 0);
  std::fill(recv_buffer, recv_buffer + 100, 0);
  std::copy_n(message.c_str(), strlen(message.c_str()) + 1, send_buffer);

  const auto l_neighbour = (comm_rank + comm_size - 1) % comm_size;
  const auto r_neighbour = (comm_rank + 1) % comm_size;

  for (auto i = 0; i < comm_size - 1; ++i) {

    if (!use_sendrecv) {
      MPI_Request send_req;
      MPI_Isend(send_buffer, 100, MPI_CHAR, r_neighbour, 27, MPI_COMM_WORLD,
                &send_req);

      MPI_Request recv_req;
      MPI_Irecv(recv_buffer, 100, MPI_CHAR, l_neighbour, 27, MPI_COMM_WORLD,
                &recv_req);

      MPI_Status send_status;
      MPI_Wait(&send_req, &send_status);

      MPI_Status recv_status;
      MPI_Wait(&recv_req, &recv_status);
    } else {
      MPI_Status sendrecv_status;
      MPI_Sendrecv(send_buffer, 100, MPI_CHAR, r_neighbour, 27, recv_buffer,
                   100, MPI_CHAR, l_neighbour, 27, MPI_COMM_WORLD,
                   &sendrecv_status);
    }

    const std::string new_message = std::string(recv_buffer) +
                                    std::string(", ") +
                                    std::to_string(comm_rank);
    std::copy_n(new_message.c_str(), strlen(new_message.c_str()) + 1,
                send_buffer);

    if (!print_only_once_per || (print_only_once_per && i == comm_size - 2)) {
      printf("%s\n", new_message.c_str());
    }
  }

  MPI_Finalize();
}
