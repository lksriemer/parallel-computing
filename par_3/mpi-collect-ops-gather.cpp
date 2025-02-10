#include <algorithm>
#include <cstring>
#include <mpi.h>
#include <random>
#include <stdio.h>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::random_device dev;
  std::seed_seq seed{static_cast<int>(dev()), rank};
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(0, 1);

  double num = dist(rng);

  int root = 0;

  auto rfield = std::vector<double>();

  if (rank == root) {
    rfield.resize(mpi_size);
  }

  MPI_Gather(&num, 1, MPI_DOUBLE, rfield.data(), 1, MPI_DOUBLE, root,
             MPI_COMM_WORLD);

  if (rank == root) {
    for (auto i = 0; i < mpi_size; ++i) {
      std::cout << rfield.at(i) << " ";
    }
    std::cout << std::endl;
  }

  MPI_Finalize();
}
