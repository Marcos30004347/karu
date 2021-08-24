#pragma once

#include "algebra/compute/Compute.hpp"

extern karu::algebra::compute::Program* blur_program; 
extern karu::algebra::compute::Kernel*  blur_kernel;

extern karu::algebra::compute::Program* resize_program; 
extern karu::algebra::compute::Kernel*  resize_kernel;

void create_sift_kernels();
void destroy_sift_kernels();