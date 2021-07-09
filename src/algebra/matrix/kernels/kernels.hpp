#pragma once

#include "algebra/core/compute/Compute.hpp"

extern karu::algebra::compute::Program* bsMV_program;
extern karu::algebra::compute::Kernel* bsMV_kernel;

void create_sparse_matrix_kernels();
void destroy_sparse_matrix_kernels();
