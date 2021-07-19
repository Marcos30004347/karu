#pragma once

#include "algebra/core/compute/Compute.hpp"

extern karu::algebra::compute::Program* radix_sort_program;
extern karu::algebra::compute::Program* array_scan_program;
extern karu::algebra::compute::Program* array_reduce_program;

extern karu::algebra::compute::Kernel*  radix_sort_int_to_int_prefix_sum_kernel;
extern karu::algebra::compute::Kernel*  radix_move_int_to_int_elements_kernel;
extern karu::algebra::compute::Kernel*  block_scan_kernel;
extern karu::algebra::compute::Kernel*  scan_kernel;
extern karu::algebra::compute::Kernel*  add_sums_kernel;
extern karu::algebra::compute::Kernel*  reduce_kernel;
extern karu::algebra::compute::Kernel*  radix_parallel_order_checking_kernel;

void create_algebra_lib_kernels();
void destroy_algebra_lib_kernels();
