#include <assert.h>
#include <iostream>

#include "algebra/core/types.hpp"
#include "algebra/compute/Compute.hpp"

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
using namespace karu::algebra::compute;


#define VECTOR_SIZE 1024

int main()
{
	Context::initContext();

	Program prog = Program(test_kernel_src);

	f32 alpha = 2.0;

	f32* A = new f32[VECTOR_SIZE];
	f32* B = new f32[VECTOR_SIZE];

	for(i32 i = 0; i < VECTOR_SIZE; i++)
	{
		A[i] = i;
		B[i] = i;
	}

	Buffer A_Buffer = Buffer(A, sizeof(f32)*VECTOR_SIZE, Buffer::READ_ONLY, false);
	Buffer B_Buffer = Buffer(B, sizeof(f32)*VECTOR_SIZE, Buffer::READ_ONLY, false);

	Buffer C_Buffer = Buffer(sizeof(f32)*VECTOR_SIZE, Buffer::WRITE_ONLY, Buffer::MEM_GPU);

	Kernel kernel = Kernel(&prog, "test_kernel");

	kernel.setKernelArgument(0, sizeof(f32), &alpha);
	kernel.setKernelArgument(1, BUFFER_ARG_SIZE, A_Buffer.upload());
	kernel.setKernelArgument(2, BUFFER_ARG_SIZE, B_Buffer.upload());
	kernel.setKernelArgument(3, BUFFER_ARG_SIZE, C_Buffer.upload());

	kernel.enqueue({ VECTOR_SIZE }, { 64 }, std::vector<Event>(), nullptr);

	// Transfer the data from the C Buffer in the Compute Unit to the Logic Unit
	C_Buffer.toLogicUnit();

	f32* data = (f32*)C_Buffer.download();

	for(i32 i = 0; i < VECTOR_SIZE; i++)
	{
		assert(data[i] == alpha * A[i] + B[i]);
	}

	Context::stopContext();
	
	delete A;
	delete B;
	
	return 0;
}
