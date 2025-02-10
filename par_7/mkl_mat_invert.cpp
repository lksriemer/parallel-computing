#include <mkl.h>

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>
#include <chrono>
#include <atomic>
#include <memory>
#include <random>

int main(int argc, char const *argv[]) {

  const int n = argc >= 2 ? std::stoi(argv[1]) : 400;
  
  std::random_device dev;
  std::seed_seq seed{static_cast<int>(dev())};
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(-100, 100);

  std::vector<double> matrix(n * n);
  std::vector<int> piv(n);
  std::vector<double> workspace(64 * n);
  const int workspace_length = 64 * n;
  std::generate(matrix.begin(), matrix.end(), [&] { return dist(rng); });

  // Warmup run
  int info;
  dgetrf(&n, &n, matrix.data(), &n, piv.data(), &info);
  dgetri(&n, matrix.data(), &n, piv.data(), workspace.data(), &workspace_length, &info);
  if(info != 0){
    std::cout << "ERROR" << std::endl;
    return 1;
  }

  const auto min_iters = 100;
  const auto min_seconds = 5;
  auto total_runtime = 0.;
  auto i = 0;
  for(; i < min_iters || (i >= min_iters && total_runtime < min_seconds * 1000000); ++i){
    const auto start = std::chrono::system_clock::now();


    std::atomic_thread_fence(std::memory_order::memory_order_seq_cst);

    // Invert
    int info;
    dgetrf(&n, &n, matrix.data(), &n, piv.data(), &info);
    dgetri(&n, matrix.data(), &n, piv.data(), workspace.data(), &workspace_length, &info);
    if(info != 0){
      std::cout << "ERROR" << std::endl;
      return 1;
    }
    // Does repeatedly (re)inverting eventually lead to singular ... ?

    std::atomic_thread_fence(std::memory_order::memory_order_seq_cst);


    const auto end = std::chrono::system_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    total_runtime += elapsed.count();
  } 

  const auto avg_runtime = total_runtime / i;

  std::cout << "avg runtime: " << avg_runtime << " at size " << n << std::endl;

  return 0;
}
