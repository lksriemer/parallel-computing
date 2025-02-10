#include <algorithm>
#include <cstring>
#include <mpi.h>
#include <random>
#include <stdio.h>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  bool do_allreduce = true;

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::random_device dev;
  std::seed_seq seed{static_cast<int>(dev()), rank};
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(0, 1);

  int root = 0;
  int pdim = 128;
  
  double local_norm = .0;
  double global_norm = .0;
  double global_norm_2 = .0;

  auto rfield = std::vector<double>();
  rfield.resize(pdim);
  std::generate(rfield.begin(), rfield.end(), [&] { return dist(rng); });

  local_norm = std::inner_product(rfield.begin(), rfield.end(), rfield.begin(), .0);

  if(do_allreduce){
    MPI_Reduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    MPI_Allreduce(&local_norm, &global_norm_2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if(root == rank && global_norm != global_norm_2){
        std::cout << "NOT SAME";
        std::cout << std::endl;
    }
  }else{
    MPI_Reduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
  }
  
  
  if (rank == root) {
    std::cout << std::sqrt(global_norm) << " ";
    std::cout << std::endl;
  }

  MPI_Finalize();
}
