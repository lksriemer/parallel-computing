#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include "args.hxx"


void transpose(const float *const input, float *const output, int n){
  for(auto i = 0; i < n; ++i){
    for(auto j = 0; j < n; ++j){
      output[i * n + j] = input[j * n + i];
    }
  }
}

void gemv_row_mpi(const float *const A_l, const float *const input, float *const output){

}

void gemv_col_comm_mpi(const float *const A_l, const float *const input, float *const output){
  
}

void gemv_col_noncomm_mpi(const float *const A_l, const float *const input, float *const output){
  
}

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  int root = 0;

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  args::ArgumentParser parser(
      "Program to solve exercise 5 of the parallel computing course");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::CompletionFlag completion(parser, {"complete"});

  args::ValueFlag<int> _n(parser, "nt", "The dimension",
                           {"n"});

  args::Flag _t_flag(parser, "t_flag",
                          "Use row (true) or col (false) partitioning",
                          {"t_flag"});
  // args::Flag _output_each(parser, "output_each",
  //                         "Output the comma-separated list at each step",
  //                         {"output_each"});
  // args::Flag _bench(parser, "bench", "Benchmark the execution time", {"bench"});

  try {
    parser.ParseCLI(argc, argv);
  } catch (const args::Completion &e) {
    std::cout << e.what();
    return 0;
  } catch (const args::Help &) {
    std::cout << parser;
    return 0;
  } catch (const args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  const auto n = args::get(_n);
  const auto t_flag = args::get(_t_flag);
  const auto bench = true;

  const auto block_size = n / mpi_size;

  if(block_size % mpi_size != 0){
    if(rank == root){
      std::cout << "mpi_size doesnt divide n" << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  const auto A_initial = std::vector<float>(n * n);
  std::fill(A_initial.begin(), A_initial.end(), 1.);
  auto A_t = std::vector<float>(n * n);
  transpose(A_initial.data(), A_t.data(), n);
  const auto A = t_flag ? A_initial : A_t;

  const auto A_l = A.data() + rank * block_size * n;

  if (bench) {
    const auto bench_runs = 10;
    for (auto bench_i = 0; bench_i < bench_runs; ++bench_i) {

      const auto start = std::chrono::system_clock::now();

      // TODO: bench here

      // Doesnt start from the same init values, but that doesnt matter for runtime
      const auto end = std::chrono::system_clock::now();
      const auto elapsed =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start);
      std::cout << elapsed.count() << ",";
    }
    std::cout << std::endl;
  }

  if (output_end) {
    for (auto i = 0; i < n_x_l; ++i) {
      std::cout << T_new->at(i) << ",";
    }
    std::cout << std::endl;
  }

  MPI_Finalize();
  return 0;
}
