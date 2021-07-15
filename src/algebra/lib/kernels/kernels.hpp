#pragma once

#include "algebra/core/compute/Compute.hpp"

extern karu::algebra::compute::Program* radix_sort_program;

extern karu::algebra::compute::Kernel*  radix_sort_int_to_int_prefix_sum_kernel;
extern karu::algebra::compute::Kernel*  radix_move_int_to_int_elements_kernel;
extern karu::algebra::compute::Kernel*  radix_scan_kernel;

void create_algebra_lib_kernels();
void destroy_algebra_lib_kernels();
