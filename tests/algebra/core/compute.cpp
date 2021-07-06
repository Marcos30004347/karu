#include <assert.h>
#include <iostream>

#include "algebra/core/compute/ComputeContext.hpp"
#include "algebra/core/compute/ComputeStorage.hpp"
#include "algebra/core/compute/ComputeKernel.hpp"
#include "algebra/core/compute/ComputeProgram.hpp"
#include "algebra/core/types.hpp"

const char* test_kernel_src =
"__kernel void test_kernel(float alpha, __global float *A, __global float *B, __global float *C) { 	\n"
"    int index = get_global_id(0);          																												\n"
"    C[index] = alpha * A[index] + B[index];																												\n"
"}                                          																												\n";

using namespace karu;
using namespace algebra;

#define VECTOR_SIZE 1024

int main()
{
	compute::Context::initContext();

	compute::Program prog = compute::Program(test_kernel_src);

	f32 alpha = 1.0;

	f32* A = new f32[VECTOR_SIZE];
	f32* B = new f32[VECTOR_SIZE];
	f32* C = new f32[VECTOR_SIZE];

	for(i32 i = 0; i < VECTOR_SIZE; i++)
	{
		A[i] = i;
		B[i] = i;
		C[i] = 0;
	}

	compute::Storage A_storage = compute::Storage((i8*)A, sizeof(f32)*VECTOR_SIZE, compute::Storage::READ_ONLY);
	compute::Storage B_storage = compute::Storage((i8*)B, sizeof(f32)*VECTOR_SIZE, compute::Storage::READ_ONLY);
	compute::Storage C_storage = compute::Storage(sizeof(f32)*VECTOR_SIZE, compute::Storage::WRITE_ONLY, compute::Storage::MEM_GPU);

	A_storage.toComputeUnit();
	B_storage.toComputeUnit();
	C_storage.toComputeUnit();

	compute::Kernel kernel = compute::Kernel(&prog, "test_kernel");
	
	kernel.setKernelArgument(0, sizeof(f32), &alpha);
	kernel.setKernelArgument(1, sizeof(cl_mem), A_storage.ref());
	kernel.setKernelArgument(2, sizeof(cl_mem), B_storage.ref());
	kernel.setKernelArgument(3, sizeof(cl_mem), C_storage.ref());

	kernel.enqueue({VECTOR_SIZE}, {64}, std::vector<compute::Event>(), nullptr);

	C_storage.toLogicUnit();

	f32* data = (f32*)C_storage.data();

	for(i32 i = 0; i < VECTOR_SIZE; i++)
	{
		std::cout << data[i] << std::endl;
	}
}
