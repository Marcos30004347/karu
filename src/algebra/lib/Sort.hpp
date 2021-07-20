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

void prefix_sum(
	Buffer* keys,
	Buffer* vals,
	int i,
	int size,
	int bs,
	Buffer* block_sum,
	Buffer* keys_shuffle,
	Buffer* vals_shuffle,
	int rounded
)
{

	// std::cout << "asdas\n";
	// int* b = (int*)keys->download();
	// for(int i=0; i<size; i++)
	// {
	// 	std::cout << b[i] << " ";
	// }
	// std::cout << "\n";
	// keys->upload();

	radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(0, keys);
	radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(1, vals);
	radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(2, &i);
	radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(3, &size);
	radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(4, sizeof(int) * roundToPowerOfTwo(bs) * 4, nullptr);
	radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(5, sizeof(int) * bs, nullptr);
	radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(6, sizeof(int) * bs, nullptr);
	radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(7, block_sum);
	radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(8, keys_shuffle);
	radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(9, vals_shuffle);
	radix_sort_int_to_int_prefix_sum_kernel->enqueue({ (unsigned long)rounded }, { (unsigned long)bs });

	// std::cout << "vals\n";
	// int* b = (int*)vals_shuffle->download();
	// for(int i=0; i<rounded; i++)
	// {
	// 	std::cout << b[i] << " ";
	// }
	// std::cout << "\n";
	// vals_shuffle->upload();

	std::cout << "keys\n";
	int* a = (int*)keys_shuffle->download();
	for(int i=0; i<rounded; i++)
	{
		std::cout << a[i] << " ";
	}
	std::cout << "\n";
	keys_shuffle->upload();


}


void move_elements(
	Buffer* keys_shuffle,
	Buffer* vals_shuffle,
	int i,
	int bs,
	Buffer* keys,
	Buffer* vals,
	Buffer* block_sum,
	Buffer* prefi_sum,
	int size,
	int rounded
)
{
	radix_move_int_to_int_elements_kernel->setKernelArgument(0, keys_shuffle);
	radix_move_int_to_int_elements_kernel->setKernelArgument(1, vals_shuffle);
	radix_move_int_to_int_elements_kernel->setKernelArgument(2, block_sum);
	radix_move_int_to_int_elements_kernel->setKernelArgument(3, prefi_sum);
	radix_move_int_to_int_elements_kernel->setKernelArgument(4, keys);
	radix_move_int_to_int_elements_kernel->setKernelArgument(5, vals);
	radix_move_int_to_int_elements_kernel->setKernelArgument(6, &i);
	radix_move_int_to_int_elements_kernel->setKernelArgument(7, &size);
	radix_move_int_to_int_elements_kernel->enqueue({ (size_t)rounded }, { (size_t)bs });
}


void sort(Buffer* keys, Buffer* vals, int size, bool do_order_check = false)
{
	int bs = 256;

	int rounded = ceil(((float)size)/bs) * bs;
	int num_blocks = rounded / bs;

	algebra::compute::Buffer block_sum = algebra::compute::Buffer(4*rounded/bs*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	algebra::compute::Buffer key_shuff = algebra::compute::Buffer(rounded*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	algebra::compute::Buffer val_shuff = algebra::compute::Buffer(rounded*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	algebra::compute::Buffer prefi_sum = algebra::compute::Buffer((4*rounded/bs)*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	algebra::compute::Buffer intermediate = algebra::compute::Buffer(num_blocks*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);

	for(int i=0; i<32; i+=2)
	{
		if(do_order_check && isInOrder(keys, rounded, num_blocks)) break;
		prefix_sum(keys, vals, i, size, bs, &block_sum, &key_shuff, &val_shuff, rounded);
		scan(&block_sum, &prefi_sum, size);
		move_elements(&key_shuff, &val_shuff, i, bs, keys, vals, &block_sum, &prefi_sum, size, rounded);
	}
}

};
