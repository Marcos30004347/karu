#pragma once

#include "algebra/compute/Compute.hpp"

extern karu::algebra::compute::Program* blur_program; 
extern karu::algebra::compute::Kernel*  blur_kernel;

void create_blur_kernels();
void destroy_blur_kernels();