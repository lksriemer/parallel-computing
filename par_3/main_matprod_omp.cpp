#include <algorithm>
#include <iostream>
#include <memory>
#include <stdio.h>
#include <vector>
#include <chrono>

#include <omp.h>

// matmuls a * b = c, nxm * mxk = nxk
void gemm_omp(std::shared_ptr<const std::vector<double>> a,
              std::shared_ptr<const std::vector<double>> b, 
              int n, int m, int k,
              std::shared_ptr<std::vector<double>> c) {
#pragma omp parallel for
  for (auto row = 0; row < n; ++row) {
    for (auto col = 0; col < k; ++col) {
      auto sum = .0;
      for (auto i = 0; i < m; ++i) {
        sum += a->at(row * m + i) * b->at(i * k + col);
      }
      c->at(row * k + col) = sum;
    }
  }
}

int main_par_3(int argc, char const *argv[]){
  const auto a = std::make_shared<std::vector<double>>([] {
    std::vector<double> v(16);
    auto n = 1.;
    std::generate(v.begin(), v.end(), [&n] { return n++; });
    return v;
  }());
  const auto b = std::make_shared<std::vector<double>>([] {
    std::vector<double> v(16);
    auto n = 16.;
    std::generate(v.begin(), v.end(), [&n] { return n--; });
    return v;
  }());
  auto c = std::make_shared<std::vector<double>>(16);

  gemm_omp(a, b, 4, 4, 4, c);

  for(auto i = 0; i < c->size(); ++i){
    std::cout << c->at(i) << " ";
  }
  std::cout << std::endl;

  return 0;
}

int main_par_7(int argc, char const *argv[]){
  const std::size_t n = argc >= 2 ? std::stoi(argv[1]) : 400;

  const auto matrix = std::make_shared<std::vector<double>>([n] {
    std::vector<double> v(n * n);
    auto n = 1.;
    std::generate(v.begin(), v.end(), [&n] { return n++; });
    return v;
  }());
  const auto input = std::make_shared<std::vector<double>>([n] {
    std::vector<double> v(n);
    auto n = 16.;
    std::generate(v.begin(), v.end(), [&n] { return n--; });
    return v;
  }());
  auto output = std::make_shared<std::vector<double>>(n);

  // Warmup run
  gemm_omp(matrix, input, n, n, 1, output);

  const auto min_iters = 100;
  const auto min_seconds = 5;
  auto total_runtime = 0.;
  auto i = 0;
  for(; i < min_iters || (i >= min_iters && total_runtime < min_seconds * 1000000); ++i){
    const auto start = std::chrono::system_clock::now();


    std::atomic_thread_fence(std::memory_order::memory_order_seq_cst);

    // Compute y = Ax
    gemm_omp(matrix, input, n, n, 1, output);

    std::atomic_thread_fence(std::memory_order::memory_order_seq_cst);


    const auto end = std::chrono::system_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    total_runtime += elapsed.count();
  } 

  const auto avg_runtime = total_runtime / i;

  std::cout << "avg runtime: " << avg_runtime << " at size " << n << std::endl;

  return 0;
}

int main(int argc, char const *argv[]) {
  return main_par_7(argc, argv);
}