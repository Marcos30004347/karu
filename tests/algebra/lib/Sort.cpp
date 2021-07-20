#include <random>
#include <chrono>
#include <assert.h>
#include <iostream>

#include "algebra/lib/Sort.hpp"
#include "algebra/core/compute/Compute.hpp"

using namespace karu;

#define N 9

int keys_test[] = {1,2,0,3,0,1,1,0,3,3,3,2,1,2,2,0,2,0,0,2};
int vals_test[] = {1,2,0,3,0,1,1,0,3,3,3,2,1,2,2,0,2,0,0,2};

unsigned sizes[N] = {
	256,
	512,
	1024,
	10000,
	100000,
	1000000,
	10000000,
	100000000,
};

std::chrono::steady_clock::time_point begin;
std::chrono::steady_clock::time_point end;

int main()
{
	algebra::compute::Context::initContext();
	// Buffer keys_buffer = Buffer(keys_test, sizeof(int)*20, algebra::compute::Buffer::READ_WRITE, false);
	// Buffer vals_buffer = Buffer(vals_test, sizeof(int)*20, algebra::compute::Buffer::READ_WRITE, false);
	
	// keys_buffer.upload();
	// vals_buffer.upload();
	
	// sort(&keys_buffer, &vals_buffer, 20);
	// int* keys = (int*)keys_buffer.download();

	// std::cout << "\n";
	// for(int j=0; j<20; j++)
	// {
	// 	std::cout << keys[j] << " ";
	// }
	// std::cout << "\n";
	// return 0;

	for(int i=0; i<1; i++)
	{
		int* keys_data = (int*)malloc(sizeof(int) * sizes[i]);
		int* vals_data = (int*)malloc(sizeof(int) * sizes[i]);
		
		for(int j=0; j<sizes[i]; j++){
        keys_data[j] = sizes[i] - j;//(rand()%sizes[i]) + 1;
        vals_data[j] = sizes[i] - j;//(rand()%sizes[i]) + 1;
		}
	
		for(int j=0; j<sizes[i]; j++)
			std::cout << keys_data[j] << " ";
		std::cout << "\n";

		Buffer keys_buffer = Buffer(keys_data, sizeof(int)*sizes[i], algebra::compute::Buffer::READ_WRITE, false);
		Buffer vals_buffer = Buffer(vals_data, sizeof(int)*sizes[i], algebra::compute::Buffer::READ_WRITE, false);
		
		keys_buffer.upload();
		vals_buffer.upload();
		
		begin = std::chrono::steady_clock::now();
		sort(&keys_buffer, &vals_buffer, sizes[i]);
		end = std::chrono::steady_clock::now();
		std::cout << "Total Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;

		int* keys = (int*)keys_buffer.download();
		int* vals = (int*)vals_buffer.download();

		std::cout << "\n";
		for(int j=0; j<sizes[i]; j++)
			std::cout << keys[j] << " ";
		std::cout << "\n";
	}


	algebra::compute::Context::stopContext();
	return 0;
}
