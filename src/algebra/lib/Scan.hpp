#include <cmath>
#include <iostream>

#include "algebra/lib/Commom.hpp"
#include "algebra/lib/kernels/kernels.hpp"
#include "algebra/core/compute/ComputeBuffer.hpp"
#include <string.h>
#include <chrono>

using namespace karu;
using namespace algebra;
using namespace compute;

namespace karu {

void scan(int* arr, int* res, int n);

void scanRec(Buffer* inp, Buffer* out, int n);
void scanBlock(Buffer* inp, Buffer* out, int n);
void scanBlocks(Buffer* inp, Buffer* out, int n);

void scanBlock(Buffer* inp, Buffer* out, int n)
{
	u32 power_of_two = roundToPowerOfTwo(n);

	scan_kernel->setKernelArgument(0, BUFFER_ARG_SIZE, inp->upload());
	scan_kernel->setKernelArgument(1, BUFFER_ARG_SIZE, out->upload());
	scan_kernel->setKernelArgument(2, 256*sizeof(int), nullptr);
	scan_kernel->setKernelArgument(3, sizeof(int), &n);
	scan_kernel->setKernelArgument(4, sizeof(int), &power_of_two);

	scan_kernel->enqueue({ power_of_two }, { power_of_two });
}

void scanBlocks(
	Buffer* inp,
	Buffer* out,
	int n
){
	int block_size = 256;

	int rounded 			= ceil(((float)n)/block_size) * block_size;
	int num_blocks 		= rounded / block_size;
	int power_of_two 	= roundToPowerOfTwo(rounded);

	Buffer sum_inp = Buffer(num_blocks*sizeof(int), Buffer::READ_WRITE, Buffer::MEM_GPU);
	Buffer sum_out = Buffer(num_blocks*sizeof(int), Buffer::READ_WRITE, Buffer::MEM_GPU);

	block_scan_kernel->setKernelArgument(0, inp);
	block_scan_kernel->setKernelArgument(1, out);
	block_scan_kernel->setKernelArgument(2, sum_inp);
	block_scan_kernel->setKernelArgument(3, block_size*sizeof(int), nullptr);
	block_scan_kernel->setKernelArgument(4, sizeof(int), &n);
	// Scan individual blocks	
	block_scan_kernel->enqueue({ (u32)rounded }, { (u32)block_size });

	// Scan the sum blocks	
	scanRec(&sum_inp, &sum_out, num_blocks);

	add_sums_kernel->setKernelArgument(0, BUFFER_ARG_SIZE, out->upload());
	add_sums_kernel->setKernelArgument(1, BUFFER_ARG_SIZE, out->upload());
	add_sums_kernel->setKernelArgument(2, BUFFER_ARG_SIZE, sum_out.upload());
	
	//Add each element of the scaned sum block to the elements of its respective block	
	add_sums_kernel->enqueue({ (u32)rounded }, { (u32)block_size });
}

void scanRec(
	Buffer* inp,
	Buffer* out,
	int n
)
{
	int block_size = 256;

	int rounded 			= ceil(((float)n)/block_size) * block_size;
	int num_blocks 		= rounded / block_size;

	if(num_blocks > 1)
		return scanBlocks(inp,out,n);

	scanBlock(inp, out, n);
}

void scan(int* arr, int* res, int n)
{
	Buffer inp = Buffer(arr, n*sizeof(int), Buffer::READ_WRITE, false);
	Buffer out = Buffer(n*sizeof(int), Buffer::READ_WRITE, Buffer::MEM_GPU);
	
	inp.toComputeUnit();
	std::chrono::steady_clock::time_point begin;
	std::chrono::steady_clock::time_point end;

	begin = std::chrono::steady_clock::now();
	scanRec(&inp, &out, n);
	end = std::chrono::steady_clock::now();

	std::cout << "GPU Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;

	// int* ret = (int*)out.download();
	// memcpy(res, ret, n*sizeof(int));
}

}
