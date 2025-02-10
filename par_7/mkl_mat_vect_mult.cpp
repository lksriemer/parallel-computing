#include <mkl.h>

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>
#include <chrono>
#include <atomic>
#include <memory>

int main(int argc, char const *argv[]) {

  const std::size_t n = argc >= 2 ? std::stoi(argv[1]) : 400;

  std::vector<double> input(n);
  std::iota(input.begin(), input.end(), 0);
  std::vector<double> output(n);
  std::fill(output.begin(), output.end(), 0);
  std::vector<double> matrix(n * n);
  std::iota(matrix.begin(), matrix.end(), 0);

  // Warmup run
  const auto alpha = 1.;
  const auto beta = 0.;
  cblas_dgemv(CBLAS_LAYOUT::CblasRowMajor, CblasNoTrans, n, n, alpha, matrix.data(), n, input.data(), 1, beta, output.data(), 1);

  const auto min_iters = 100;
  const auto min_seconds = 5;
  auto total_runtime = 0.;
  auto i = 0;
  for(; i < min_iters || (i >= min_iters && total_runtime < min_seconds * 1000000); ++i){
    const auto start = std::chrono::system_clock::now();


    std::atomic_thread_fence(std::memory_order::memory_order_seq_cst);

    // Compute y = Ax
    cblas_dgemv(CBLAS_LAYOUT::CblasRowMajor, CblasNoTrans, n, n, alpha, matrix.data(), n, input.data(), 1, beta, output.data(), 1);

    std::atomic_thread_fence(std::memory_order::memory_order_seq_cst);


    const auto end = std::chrono::system_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    total_runtime += elapsed.count();
  } 

  const auto avg_runtime = total_runtime / i;

  std::cout << "avg runtime: " << avg_runtime << " at size " << n << std::endl;

  return 0;
}
