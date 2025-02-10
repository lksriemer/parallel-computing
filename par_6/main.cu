#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include <cub/cub.cuh>

#include "cublas_v2.h"
#include <cuda_runtime.h>

constexpr int block_size = 256;

__global__ void vectorScalarProduct(float *vector, float scalar,
                                    std::size_t n) {
  const auto gid = threadIdx.x + blockIdx.x * blockDim.x;

  if (gid >= n) {
    return;
  }

  vector[gid] *= scalar;
}

__global__ void gemv(float *m, float *input, float *output, std::size_t n) {
  using BlockReduce = cub::BlockReduce<float, block_size>;
  __shared__ typename BlockReduce::TempStorage temp_storage;

  const auto gid = threadIdx.x + blockIdx.x * blockDim.x;
  const auto tid = threadIdx.x;

  const auto block_start = gid / block_size * block_size;
  const auto block_end = block_start + block_size - 1;

  if (gid >= n * n) {
    return;
  }

  const float v = m[gid] * input[gid % n];

  const auto row = gid / n;

  const auto start_row = block_start / n;
  const auto end_row = block_end / n;

  for (auto cur_row = start_row; cur_row <= end_row; ++cur_row) {
    const float sum = BlockReduce(temp_storage).Sum(row == cur_row ? v : .0);
    if (tid == 0) {
      atomicAdd(output + cur_row, sum);
    }
  }
}

__global__ void gemv_row_blocked(float *mat, float *input, float *output,
                                 std::size_t n, std::size_t m) {
  using BlockReduce = cub::BlockReduce<float, block_size>;
  __shared__ typename BlockReduce::TempStorage temp_storage;

  const auto gid = threadIdx.x + blockIdx.x * blockDim.x;
  const auto tid = threadIdx.x;

  if (gid >= block_size * m) {
    return;
  }

  const auto row = gid / block_size;

  float sum = .0;
  for (auto i = tid; i < n; i += block_size) {
    sum += mat[row * n + i] * input[i];
  }

  // Return the warp-wide sum to lane0
  const float full_sum = BlockReduce(temp_storage).Sum(sum);
  if (tid == 0) {
    output[row] = full_sum;
  }
}

__global__ void gemv_row_blocked_d(double *mat, double *input, double *output,
                                   std::size_t n, std::size_t m) {
  using BlockReduce = cub::BlockReduce<double, block_size>;
  __shared__ typename BlockReduce::TempStorage temp_storage;

  const auto gid = threadIdx.x + blockIdx.x * blockDim.x;
  const auto tid = threadIdx.x;

  if (gid >= block_size * m) {
    return;
  }

  const auto row = gid / block_size;

  double sum = .0;
  for (auto i = tid; i < n; i += block_size) {
    sum += mat[row * n + i] * input[i];
  }

  // Return the warp-wide sum to lane0
  const double full_sum = BlockReduce(temp_storage).Sum(sum);
  if (tid == 0) {
    output[row] = full_sum;
  }
}

void gen_rand_vec(float *v, int n) {
  std::uniform_real_distribution<float> distribution(-1000., 1000.);
  std::mt19937 engine;
  auto generator = std::bind(distribution, engine);
  std::generate_n(v, n, generator);
}

void gemv_cpu(float *mat, float *input, float *output, std::size_t n,
              std::size_t m) {
  for (auto row = 0; row < m; ++row) {
    for (auto col = 0; col < n; ++col) {
      output[row] += mat[row * n + col] * input[col];
    }
  }
}

int main_par_6(int argc, char *argv[]) {
  const std::size_t n = argc >= 2 ? std::stoi(argv[1]) : 10;
  const std::size_t m = argc >= 3 ? std::stoi(argv[2]) : n;

  cublasHandle_t handle;
  cublasCreate(&handle);

  std::vector<float> h_input(n);
  std::iota(h_input.begin(), h_input.end(), 0);
  std::vector<float> h_output(n);
  std::fill(h_output.begin(), h_output.end(), 0);
  std::vector<float> h_matrix(n * n);
  std::iota(h_matrix.begin(), h_matrix.end(), 0);

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  float *d_input;
  cudaMalloc(&d_input, n * sizeof(float));
  cudaMemcpy(d_input, h_input.data(), n * sizeof(float),
             cudaMemcpyHostToDevice);
  float *d_output;
  cudaMalloc(&d_output, n * sizeof(float));
  cudaMemcpy(d_output, h_output.data(), n * sizeof(float),
             cudaMemcpyHostToDevice);
  float *d_matrix;
  cudaMalloc(&d_matrix, n * n * sizeof(float));
  cudaMemcpy(d_matrix, h_matrix.data(), n * n * sizeof(float),
             cudaMemcpyHostToDevice);

  int blockSize = block_size;
  int gridSize = (m * block_size + blockSize - 1) / blockSize;

  float one = 1.;
  float zero = 0.;

  cudaEventRecord(start);

  //   cublasSgemv(handle, CUBLAS_OP_T,
  //     m, n,
  //     &one,
  //     d_matrix, m,
  //     d_input, 1,
  //     &zero,
  //     d_output, 1);
  gemv_row_blocked<<<gridSize, blockSize>>>(d_matrix, d_input, d_output, n, m);

  cudaEventRecord(stop);

  cudaMemcpy(h_output.data(), d_output, n * sizeof(float),
             cudaMemcpyDeviceToHost);

  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);

  //   for (int i = 0; i < n; ++i) {
  //     std::cout << h_output[i] << " ";
  //   }
  //   std::cout << std::endl;

  const float ops = n * n + 2 * n;
  std::cout << "took " << milliseconds << "ms which is "
            << (ops * sizeof(float)) / milliseconds * 1000. / 1000000000.
            << "GByte/s" << std::endl;

  return 0;
}

int main_par_7(int argc, char *argv[]) {
  const std::size_t n = argc >= 2 ? std::stoi(argv[1]) : 10;
  const std::size_t m = argc >= 3 ? std::stoi(argv[2]) : n;

  cublasHandle_t handle;
  cublasCreate(&handle);

  std::vector<double> h_input(n);
  std::iota(h_input.begin(), h_input.end(), 0);
  std::vector<double> h_output(n);
  std::fill(h_output.begin(), h_output.end(), 0);
  std::vector<double> h_matrix(n * n);
  std::iota(h_matrix.begin(), h_matrix.end(), 0);

  double *d_input;
  cudaMalloc(&d_input, n * sizeof(double));
  cudaMemcpy(d_input, h_input.data(), n * sizeof(double),
             cudaMemcpyHostToDevice);
  double *d_output;
  cudaMalloc(&d_output, n * sizeof(double));
  cudaMemcpy(d_output, h_output.data(), n * sizeof(double),
             cudaMemcpyHostToDevice);
  double *d_matrix;
  cudaMalloc(&d_matrix, n * n * sizeof(double));
  cudaMemcpy(d_matrix, h_matrix.data(), n * n * sizeof(double),
             cudaMemcpyHostToDevice);

  const auto min_iters = 100;
  const auto min_seconds = 5;
  auto total_runtime = 0.;
  auto i = 0;
  for (; i < min_iters ||
         (i >= min_iters && total_runtime < min_seconds * 1000000);
       ++i) {

    // Compute y = Ax
    int blockSize = block_size;
    int gridSize = (m * block_size + blockSize - 1) / blockSize;

    double one = 1.;
    double zero = 0.;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);

      cublasDgemv(handle, CUBLAS_OP_T,
        m, n,
        &one,
        d_matrix, m,
        d_input, 1,
        &zero,
        d_output, 1);
    // gemv_row_blocked_d<<<gridSize, blockSize>>>(d_matrix, d_input, d_output, n,
    //                                           m);

    cudaEventRecord(stop);

    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    const auto elapsed = (int)(milliseconds * 1000);
    total_runtime += elapsed;
  }

  const auto avg_runtime = total_runtime / i;

  std::cout << "avg runtime: " << avg_runtime << " at size " << n << std::endl;

  return 0;
}

int main(int argc, char *argv[]) { return main_par_7(argc, argv); }