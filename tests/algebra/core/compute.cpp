#include <assert.h>
#include <iostream>

#include "algebra/core/types.hpp"
#include "algebra/core/compute/Compute.hpp"

const char* test_kernel_src =
"__kernel void test_kernel(\n"
"	float alpha,\n"
"	__global float *A,\n"
"	__global float *B,\n"
"	__global float *C\n"
") {\n"
"    int index = get_global_id(0);\n"
"    C[index] = alpha * A[index] + B[index];\n"
"}\n";

using namespace karu;
using namespace algebra;

#define VECTOR_SIZE 1024

int main()
{
	compute::Context::initContext();

	compute::Program prog = compute::Program(test_kernel_src);

	f32 alpha = 2.0;

	f32* A = new f32[VECTOR_SIZE];
	f32* B = new f32[VECTOR_SIZE];

	for(i32 i = 0; i < VECTOR_SIZE; i++)
	{
		A[i] = i;
		B[i] = i;
	}

	compute::Buffer A_Buffer = compute::Buffer(A, sizeof(f32)*VECTOR_SIZE, compute::Buffer::READ_ONLY, false);
	compute::Buffer B_Buffer = compute::Buffer(B, sizeof(f32)*VECTOR_SIZE, compute::Buffer::READ_ONLY, false);

	compute::Buffer C_Buffer = compute::Buffer(sizeof(f32)*VECTOR_SIZE, compute::Buffer::WRITE_ONLY, compute::Buffer::MEM_GPU);

	compute::Kernel kernel = compute::Kernel(&prog, "test_kernel");

	kernel.setKernelArgument(0, sizeof(f32), &alpha);
	kernel.setKernelArgument(1, BUFFER_ARG_SIZE, A_Buffer.upload());
	kernel.setKernelArgument(2, BUFFER_ARG_SIZE, B_Buffer.upload());
	kernel.setKernelArgument(3, BUFFER_ARG_SIZE, C_Buffer.upload());

	kernel.enqueue({ VECTOR_SIZE }, { 64 }, std::vector<compute::Event>(), nullptr);

	// Transfer the data from the C Buffer in the Compute Unit to the Logic Unit
	C_Buffer.toLogicUnit();

	f32* data = (f32*)C_Buffer.download();

	for(i32 i = 0; i < VECTOR_SIZE; i++)
	{
		assert(data[i] == alpha * A[i] + B[i]);
	}

	compute::Context::stopContext();
	
	delete A;
	delete B;
	
	return 0;
}
