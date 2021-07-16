#include <cmath>
#include <iostream>

#include "algebra/lib/Commom.hpp"
#include "algebra/lib/kernels/kernels.hpp"
#include "algebra/core/compute/ComputeBuffer.hpp"
#include <string.h>

using namespace karu;
using namespace algebra;
using namespace compute;

namespace karu {

// TODO: refactoring, work sizes may not be correct

void scan_small(Buffer* inp, Buffer* out,int n)
{
	u32 power_of_two = roundToPowerOfTwo(n);
	scan_kernel->setKernelArgument(0, BUFFER_ARG_SIZE, inp->computeUnitRef());
	scan_kernel->setKernelArgument(1, BUFFER_ARG_SIZE, out->computeUnitRef());
	scan_kernel->setKernelArgument(2, power_of_two*sizeof(int), nullptr);
	scan_kernel->setKernelArgument(3, sizeof(int), &n);
	scan_kernel->setKernelArgument(4, sizeof(int), &power_of_two);
	scan_kernel->enqueue({ power_of_two }, { power_of_two });
}

void scan(int* arr, int* res, int n);
void scan_rec(Buffer* inp, Buffer* out, int n);

void scan_blocks(
	Buffer* inp,
	Buffer* out,
	int n
){
	int block_size = 4;

	int rounded 			= ceil(((float)n)/block_size) * block_size;
	int num_blocks 		= rounded / block_size;
	int power_of_two 	= roundToPowerOfTwo(rounded);

	Buffer sum_inp = Buffer(num_blocks*sizeof(int), Buffer::READ_WRITE, Buffer::MEM_GPU);
	Buffer sum_out = Buffer(num_blocks*sizeof(int), Buffer::READ_WRITE, Buffer::MEM_GPU);

	block_scan_kernel->setKernelArgument(0, BUFFER_ARG_SIZE, inp->computeUnitRef());
	block_scan_kernel->setKernelArgument(1, BUFFER_ARG_SIZE, out->computeUnitRef());
	block_scan_kernel->setKernelArgument(2, BUFFER_ARG_SIZE, sum_inp.computeUnitRef());
	block_scan_kernel->setKernelArgument(3, power_of_two*sizeof(int), nullptr);
	block_scan_kernel->setKernelArgument(4, sizeof(int), &n);

	// Scan individual blocks	
	block_scan_kernel->enqueue({ (u32)power_of_two }, { (u32)block_size });

	// Scan the sum blocks	
	scan_rec(&sum_inp, &sum_out, num_blocks);

	add_sums_kernel->setKernelArgument(0, BUFFER_ARG_SIZE, out->computeUnitRef());
	add_sums_kernel->setKernelArgument(1, BUFFER_ARG_SIZE, out->computeUnitRef());
	add_sums_kernel->setKernelArgument(2, BUFFER_ARG_SIZE, sum_out.computeUnitRef());
	
	// Add each element of the scaned sum block to the elements of its respective block	
	add_sums_kernel->enqueue({ (u32)power_of_two }, { (u32)block_size });
}

void scan_rec(
	Buffer* inp,
	Buffer* out,
	int n
)
{
	int block_size = 4;

	int rounded 			= ceil(((float)n)/block_size) * block_size;
	int num_blocks 		= rounded / block_size;

	if(num_blocks > 1)
		return scan_blocks(inp,out,n);

	scan_small(inp, out, n);
}

void scan(int* arr, int* res, int n)
{
	std::cout << "input: ";
	for(int i=0;i<n; i++)
		std::cout << arr[i] << " ";
	std::cout << std::endl;

	Buffer inp = Buffer(arr, n*sizeof(int), Buffer::READ_WRITE, false);
	Buffer out = Buffer(n*sizeof(int), Buffer::READ_WRITE, Buffer::MEM_GPU);

	scan_rec(&inp, &out, n);

	int* ret = (int*)out.logicUnitRef();
	memcpy(res, ret, n*sizeof(int));

	std::cout << "output: ";
	for(int i=0; i<n; i++)
	{
		std::cout << res[i] << " ";
	}
	std::cout << std::endl;
}

}
