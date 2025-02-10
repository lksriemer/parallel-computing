#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include "args.hxx"


void output_t_mpi(std::shared_ptr<std::vector<double>> T, int mpi_size, int rank, int root){
  const auto n_x_l = T->size();
  auto T_full = std::vector<double>();
  if(rank == root){
    T_full.resize(n_x_l * mpi_size);
  }

  MPI_Gather(T->data(), n_x_l, MPI_DOUBLE, T_full.data(), n_x_l, MPI_DOUBLE, root, MPI_COMM_WORLD);

  if(rank == root){
    for (auto i = 0; i < n_x_l * mpi_size; ++i) {
      std::cout << T_full.at(i) << ",";
    }
    std::cout << std::endl;
  }
}


void compute_temp(std::shared_ptr<std::vector<double>> T_old, 
                  std::shared_ptr<std::vector<double>> T_new, 
                  double factor, int n_x_l, int n_t, 
                  int rank, int root, int mpi_size, 
                  bool output_each){
  for (auto t = 0; t < n_t; ++t) {

    if(output_each){
      output_t_mpi(T_old, mpi_size, rank, root)
    }
    // TODO: Start send for left and right ghost values (if applicable)
    
      MPI_Isend(T_old->data())

      // TODO: Start receive for left and right ghost values (if applicable)
    
        for (auto i = 1; i < n_x_l - 1; ++i) {
          T_new->at(i) =
              T_old->at(i) +
              factor * (T_old->at(i - 1) - 2 * T_old->at(i) + T_old->at(i + 1));
        }

        // TODO: Finish receive on left and right ghost values (if applicable)

        // TODO: Compute on border values (if applicable)

        // TODO: Finish send on left and right ghost values (if applicable)

        std::copy_n(T_new->begin(), n_x_l, T_old->begin());
  }
}

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  int root = 0;

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  args::ArgumentParser parser(
      "Program to solve exercise 4 of the parallel computing course");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::CompletionFlag completion(parser, {"complete"});

  args::ValueFlag<double> t_0(parser, "t0", "The left temperature", {"t0"});
  args::ValueFlag<double> t_1(parser, "t1", "The right temperature", {"t1"});
  args::ValueFlag<double> t_i(parser, "ti", "The initial middle temperature",
                              {"ti"});
  args::ValueFlag<double> _width(parser, "width", "The width", {"w"});
  args::ValueFlag<double> _lambda(parser, "lambda", "The thermal conductivity",
                                  {"l"});
  args::ValueFlag<double> _time_end(parser, "time_end", "The end time", {"e"});
  args::ValueFlag<int> _nx(parser, "nx", "The number of space discrete steps",
                           {"nx"});
  args::ValueFlag<int> _nt(parser, "nt", "The number of time discrete steps",
                           {"nt"});

  args::Flag _output_end(parser, "output_end",
                          "Output the comma-separated list when done",
                          {"output_end"});
  args::Flag _output_each(parser, "output_each",
                          "Output the comma-separated list at each step",
                          {"output_each"});
  args::Flag _bench(parser, "bench", "Benchmark the execution time", {"bench"});

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

  const auto output_end = args::get(_output_end);
  const auto output_each = args::get(_output_each);
  const auto bench = args::get(_bench);

  const auto left_temp = args::get(t_0);
  const auto right_temp = args::get(t_1);
  const auto init_temp = args::get(t_i);

  const auto width = args::get(_width); // s

  const auto lambda = args::get(_lambda); // thermal conductivity

  const auto t_end = args::get(_time_end);

  const auto n_t = args::get(_nt);
  const auto n_x = args::get(_nx);

  if(n_x % mpi_size != 0){
    if(rank == root){
      std::cout << "mpi_size doesnt divide n_x" << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  const auto n_x_l = n_x / mpi_size;

  const auto delta_t = t_end / (double)n_t;
  const auto delta_x = width / (double)n_x;

  const auto T_old = std::make_shared<std::vector<double>>(n_x_l);
  const auto T_new = std::make_shared<std::vector<double>>(n_x_l);
  std::fill(T_old->begin(), T_old->end(), init_temp);

  const auto l_neighbour = (rank + mpi_size - 1) % mpi_size;
  const auto r_neighbour = (rank + 1) % mpi_size;

  // Init left and right temp
  if(rank == 0){
    T_old->at(0) = left_temp;
    T_new->at(0) = left_temp;
  }
  
  if(rank == mpi_size - 1){
    T_old->at(n_x_l - 1) = right_temp;
    T_new->at(n_x_l - 1) = right_temp;
  }
  

  const auto factor = lambda * delta_t / (delta_x * delta_x);

  // Can be seen as warmup run
  compute_temp(T_old, T_new, factor, n_x_l, n_t, rank, root, mpi_size, output_each);

  if (bench) {
    const auto bench_runs = 10;
    for (auto bench_i = 0; bench_i < bench_runs; ++bench_i) {

      const auto start = std::chrono::system_clock::now();

      compute_temp(T_old, T_new, factor, n_x_l, n_t, rank, root, mpi_size, output_each);
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
