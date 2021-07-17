#include <cmath>
#include <iostream>

#include "algebra/lib/Commom.hpp"
#include "algebra/lib/kernels/kernels.hpp"
#include "algebra/core/compute/ComputeBuffer.hpp"

namespace karu {

void sort(
	int* keys,
	int* vals,
	int size,
	int bits_count
)
{
	int bs = 4;

	int closesMultile = ceil(((float)size)/bs) * bs;

	int next_multiple_of_two = roundToPowerOfTwo(size);

	int num_blocks = closesMultile / bs;

	algebra::compute::Buffer keys_io = algebra::compute::Buffer(keys, size*sizeof(int), algebra::compute::Buffer::READ_WRITE, false);
	algebra::compute::Buffer vals_io = algebra::compute::Buffer(vals, size*sizeof(int), algebra::compute::Buffer::READ_WRITE, false);

	algebra::compute::Buffer block_sum = algebra::compute::Buffer(roundToPowerOfTwo(4*closesMultile/bs)*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	algebra::compute::Buffer key_shuff = algebra::compute::Buffer(roundToPowerOfTwo(closesMultile)*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	algebra::compute::Buffer val_shuff = algebra::compute::Buffer(roundToPowerOfTwo(closesMultile)*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	algebra::compute::Buffer prefi_sum = algebra::compute::Buffer(roundToPowerOfTwo(4 * size/bs)*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);
	
	algebra::compute::Buffer intermediate = algebra::compute::Buffer(num_blocks*sizeof(int), algebra::compute::Buffer::READ_WRITE, algebra::compute::Buffer::MEM_GPU);

	for(int i=0; i<2; i+=2)
	{
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(0, BUFFER_ARG_SIZE, keys_io.upload());
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(1, BUFFER_ARG_SIZE, vals_io.upload());

		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(2, sizeof(int), &i);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(3, sizeof(int), &size);

		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(4, sizeof(int) * bs * 4, nullptr);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(5, sizeof(int) * bs, nullptr);
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(6, sizeof(int) * bs, nullptr);

		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(7, BUFFER_ARG_SIZE, block_sum.upload());
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(8, BUFFER_ARG_SIZE, key_shuff.upload());
		radix_sort_int_to_int_prefix_sum_kernel->setKernelArgument(9, BUFFER_ARG_SIZE, val_shuff.upload());

		// Dispatch kernel
		radix_sort_int_to_int_prefix_sum_kernel->enqueue({ (unsigned long)closesMultile }, { (unsigned long)bs });

		std::cout << "key: ";
		int* keys_ptr = (int*)keys_io.download();
		std::cout << std::endl;
		for(int i=0; i<size;i++)
		{
			std::cout << keys_ptr[i]<< " ";
		}
		std::cout << std::endl;

		std::cout << "local shuffle: ";
		int* key_shu = (int*)key_shuff.download();
		std::cout << std::endl;
		for(int i=0; i<size;i++)
		{
			std::cout << key_shu[i]<< " ";
		}
		std::cout << std::endl;

		// std::cout << " local shuffle: ";
		// int* local_pref_shu = (int*)val_shuff.download();
		// std::cout << std::endl;
		// for(int i=0; i<size;i++)
		// {
		// 	if(i%bs == 0)
		// 		std::cout << " ";
		// 	std::cout << local_pref_shu[i]<< " ";
		// }
		// std::cout << std::endl;
	
		std::cout << "block sum: ";
		int* block_sum_ptr = (int*)block_sum.download();
		std::cout << std::endl;
		for(int i=0; i<size;i++)
		{
			std::cout << block_sum_ptr[i]<< " ";
		}
		std::cout << std::endl;

		// int wgSize = std::min(
		// 	radix_scan_kernel->getWorkGroupSize(),
		// 	(size_t)closesMultile
		// );

		// int localSize = std::min(wgSize, closesMultile);
		
		// radix_scan_kernel->setKernelArgument(2, BUFFER_ARG_SIZE, intermediate.upload());
		// radix_scan_kernel->setKernelArgument(0, BUFFER_ARG_SIZE, prefi_sum.upload());
		// radix_scan_kernel->setKernelArgument(1, BUFFER_ARG_SIZE, block_sum.upload());
		// // radix_scan_kernel->setKernelArgument(2, BUFFER_ARG_SIZE, intermediate.upload());
		// radix_scan_kernel->setKernelArgument(2, (closesMultile)*sizeof(int), nullptr);
		// radix_scan_kernel->setKernelArgument(3, sizeof(int), &size);
		// radix_scan_kernel->setKernelArgument(4, sizeof(int), &next_multiple_of_two);
		// radix_scan_kernel->setKernelArgument(4, sizeof(int), &closesMultile);
	
		// radix_scan_kernel->enqueue({(unsigned long)((size+1)/2)}, {1});
		// radix_scan_kernel->enqueue({(unsigned long)(size)}, {1});
		
		std::cout << "prefix sum: ";
		int* pref_ptr = (int*)prefi_sum.download();
		std::cout << std::endl;
		for(int i=0; i<size;i++)
		{
			std::cout << pref_ptr[i]<< " ";
		}
		std::cout << std::endl;
		// std::cout << std::endl;

		// radix_move_int_to_int_elements_kernel->setKernelArgument(0, BUFFER_ARG_SIZE, key_shuff.upload());
		// radix_move_int_to_int_elements_kernel->setKernelArgument(1, BUFFER_ARG_SIZE, val_shuff.upload());
		// radix_move_int_to_int_elements_kernel->setKernelArgument(2, 4 * bs * sizeof(int), nullptr);
		// radix_move_int_to_int_elements_kernel->setKernelArgument(3, 4 * bs * sizeof(int), nullptr);
		// radix_move_int_to_int_elements_kernel->setKernelArgument(4, BUFFER_ARG_SIZE, block_sum.upload());
		// radix_move_int_to_int_elements_kernel->setKernelArgument(5, BUFFER_ARG_SIZE, prefi_sum.upload());
		// radix_move_int_to_int_elements_kernel->setKernelArgument(6, BUFFER_ARG_SIZE, keys_io.upload());
		// radix_move_int_to_int_elements_kernel->setKernelArgument(7, BUFFER_ARG_SIZE, vals_io.upload());
		// radix_move_int_to_int_elements_kernel->setKernelArgument(8, sizeof(int), &i);
		// radix_move_int_to_int_elements_kernel->setKernelArgument(9, sizeof(int), &closesMultile);

		// radix_move_int_to_int_elements_kernel->enqueue({ (unsigned long)closesMultile }, { (unsigned long)bs });
	}

	int* keys_ptr = (int*)keys_io.download();
	int* vals_ptr = (int*)vals_io.download();

	std::cout << std::endl;
	std::cout <<  "keys: " << std::endl;
	for(int i=0; i<size;i++)
		std::cout << keys_ptr[i]<< " ";
	std::cout << std::endl;

}

};
