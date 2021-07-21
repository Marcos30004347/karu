
#include "algebra/compute/lib/Commom.hpp"
#include "algebra/compute/Buffer.hpp"
#include "algebra/compute/lib/kernels/kernels.hpp"
#include <iostream>
#include <string.h>
#include <math.h>


namespace karu::algebra::compute {

void reduce(Buffer* inp, Buffer* out, size_t n)
{
	int block_size = 256;
	int rounded 			= ceil(((float)n)/block_size) * block_size;
	int num_blocks 		= rounded / block_size;
	int power_of_two 	= roundToPowerOfTwo(rounded);

	int offset = 0;

	Buffer* output = out;

	if(num_blocks > 1)
		output = new Buffer(num_blocks*sizeof(int), Buffer::READ_WRITE, Buffer::MEM_GPU);

	reduce_kernel->setKernelArgument(0, inp);
	reduce_kernel->setKernelArgument(1, output);
	reduce_kernel->setKernelArgument(2, block_size*sizeof(int), nullptr);
	reduce_kernel->setKernelArgument(3, offset);

	reduce_kernel->enqueue({ (u32)rounded }, { (u32)block_size });

	if(num_blocks > 1)
	{
		reduce(output, out, num_blocks);
		delete output;
	}
}

// void reduceRec(
// 	Buffer* inp,
// 	Buffer* out,
// 	int n
// )
// {
// 	int block_size = 256;

// 	int rounded 			= ceil(((float)n)/block_size) * block_size;
// 	int num_blocks 		= rounded / block_size;

// 	if(num_blocks > 1)
// 		return reduceBlocks(inp,out,n);
// }

// void reduceInt32(Buffer* inp, Buffer* out, size_t n)
// {
// 	reduceInt32Blocks(inp, out, n);
// 	// int* ret = (int*)out.download();
// 	// memcpy(res, ret, n*sizeof(int));
// }

}
