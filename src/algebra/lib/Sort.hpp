#include <cmath>
#include <iostream>

#include "algebra/lib/Scan.hpp"
#include "algebra/lib/Commom.hpp"
#include "algebra/lib/Reduce.hpp"
#include "algebra/lib/kernels/kernels.hpp"
#include "algebra/core/compute/ComputeBuffer.hpp"

namespace karu {

bool isInOrder(Buffer* input, int size, int num_blocks)
{
	Buffer output = Buffer(num_blocks*sizeof(int), Buffer::READ_WRITE, Buffer::MEM_GPU);

	radix_parallel_order_checking_kernel->setKernelArgument(0, input);
	radix_parallel_order_checking_kernel->setKernelArgument(1, &output);
	radix_parallel_order_checking_kernel->setKernelArgument(2, num_blocks*sizeof(int), nullptr);
	radix_parallel_order_checking_kernel->enqueue({ (u32)size }, { (u32)(size/num_blocks) });

	Buffer result = Buffer(sizeof(int), Buffer::READ_WRITE, Buffer::MEM_GPU);

	reduceInt32(&output, &result, num_blocks);

	return !(*(int*)result.download()) > 0;
}

void sort(Buffer* keys, Buffer* vals, int size, bool do_order_check = false)
{
	int bs = 4;
	int rounded = ceil(((float)size)/bs) * bs;
	int next_multiple_of_two = karu::roundToPowerOfTwo(size);
	int num_blocks = rounded / bs;

	algebra::compute::Buffer block_sum = algebra::compute::Buffer(roundToPowerOfTwo(4*rounded/bs)*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	algebra::compute::Buffer key_shuff = algebra::compute::Buffer(roundToPowerOfTwo(rounded)*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	algebra::compute::Buffer val_shuff = algebra::compute::Buffer(roundToPowerOfTwo(rounded)*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	algebra::compute::Buffer prefi_sum = algebra::compute::Buffer(roundToPowerOfTwo(4 * size/bs)*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	
	algebra::compute::Buffer intermediate = algebra::compute::Buffer(num_blocks*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);

	for(int i=0; i<32; i+=2)
	{
		if(do_order_check && isInOrder(keys, rounded, num_blocks)) break;
		
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(0, keys);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(1, vals);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(2, &i);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(3, &size);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(4, sizeof(int) * bs * 4, nullptr);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(5, sizeof(int) * bs, nullptr);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(6, sizeof(int) * bs, nullptr);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(7, &block_sum);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(8, &key_shuff);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(9, &val_shuff);
		radix_sort_int_to_int_prefix_sum_kernel->enqueue({ (unsigned long)rounded }, { (unsigned long)bs });

		scan(&block_sum, &prefi_sum, size);

	
		radix_move_int_to_int_elements_kernel->setKernelArgument(0, &key_shuff);
		radix_move_int_to_int_elements_kernel->setKernelArgument(1, &val_shuff);
		radix_move_int_to_int_elements_kernel->setKernelArgument(2, 4 * bs * sizeof(int), nullptr);
		radix_move_int_to_int_elements_kernel->setKernelArgument(3, 4 * bs * sizeof(int), nullptr);
		radix_move_int_to_int_elements_kernel->setKernelArgument(4, block_sum);
		radix_move_int_to_int_elements_kernel->setKernelArgument(5, prefi_sum);
		radix_move_int_to_int_elements_kernel->setKernelArgument(6, keys);
		radix_move_int_to_int_elements_kernel->setKernelArgument(7, vals);
		radix_move_int_to_int_elements_kernel->setKernelArgument(8, &i);
		radix_move_int_to_int_elements_kernel->setKernelArgument(9, &rounded);
		radix_move_int_to_int_elements_kernel->enqueue({ (size_t)rounded }, { (size_t)bs });
	}
}

};
