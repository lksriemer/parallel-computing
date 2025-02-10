#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>
#include <chrono>

#include "args.hxx"

int main(int argc, char const *argv[]) {

  args::ArgumentParser parser(
      "Program to solve exercise 2 of the parallel computing course");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::CompletionFlag completion(parser, {"complete"});

  args::ValueFlag<double> t_0(parser, "t0", "The left temperature", {"t0"});
  args::ValueFlag<double> t_1(parser, "t1", "The right temperature", {"t1"});
  args::ValueFlag<double> t_i(parser, "ti", "The initial middle temperature",
                              {"ti"});
  args::ValueFlag<double> _width(parser, "width", "The width", {"w"});
  args::ValueFlag<double> _lambda(parser, "lambda", "The thermal conductivity", {"l"});
  args::ValueFlag<double> _time_end(parser, "time_end", "The end time", {"e"});
  args::ValueFlag<int> _nx(parser, "nx", "The number of space discrete steps", {"nx"});
  args::ValueFlag<int> _nt(parser, "nt", "The number of time discrete steps", {"nt"});

  args::Flag _output_list(parser, "output_list",
                          "Output the comma-separated list when done",
                          {"output_list"});
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

  const auto output_list = args::get(_output_list);
  const auto bench = args::get(_bench);

  const auto left_temp = args::get(t_0);
  const auto right_temp = args::get(t_1);
  const auto init_temp = args::get(t_i);

  const auto width = args::get(_width); // s

  const auto lambda = args::get(_lambda); // thermal conductivity

  const auto t_end = args::get(_time_end);

  const auto n_t = args::get(_nt);
  const auto n_x = args::get(_nx);

  const auto delta_t = t_end / (double)n_t;
  const auto delta_x = width / (double)n_x;

  const auto T_old = std::make_shared<std::vector<double>>(n_x);
  const auto T_new = std::make_shared<std::vector<double>>(n_x);
  std::fill(T_old->begin(), T_old->end(), init_temp);

  // Init left and right temp
  T_old->at(0) = left_temp;
  T_old->at(n_x - 1) = right_temp;
  T_new->at(0) = left_temp;
  T_new->at(n_x - 1) = right_temp;

  const auto factor = lambda * delta_t / (delta_x * delta_x);

  // Can be seen as warmup run
  for (auto t = 0; t < n_t; ++t) {
    for (auto i = 1; i < n_x - 1; ++i) {
      T_new->at(i) =
          T_old->at(i) +
          factor * (T_old->at(i - 1) - 2 * T_old->at(i) + T_old->at(i + 1));
    }

    std::copy_n(T_new->begin(), n_x, T_old->begin());
  }

  if (bench) {
    const auto bench_runs = 10;
    for (auto bench_i = 0; bench_i < bench_runs; ++bench_i) {

      const auto start = std::chrono::system_clock::now();

      for (auto t = 0; t < n_t; ++t) {
        for (auto i = 1; i < n_x - 1; ++i) {
          T_new->at(i) =
              T_old->at(i) +
              factor * (T_old->at(i - 1) - 2 * T_old->at(i) + T_old->at(i + 1));
        }
        std::copy_n(T_new->begin(), n_x, T_old->begin());
      }
      // Doesnt start from the same init values, but that doesnt matter for runtime
      const auto end = std::chrono::system_clock::now();
      const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
      std::cout << elapsed.count() << ",";
    }
    std::cout << std::endl;
  }

  if (output_list) {
    for (auto i = 0; i < n_x; ++i) {
      std::cout << T_new->at(i) << ",";
    }
    std::cout << std::endl;
  }

  return 0;
}
