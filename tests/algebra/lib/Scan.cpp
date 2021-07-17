#include <chrono>
#include <random>
#include <assert.h>
#include <iostream>

#include "algebra/lib/Scan.hpp"
#include "algebra/core/compute/Compute.hpp"

using namespace karu;

// Number of tests to execute
#define N 9

// The array lenght of each test
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

// Naive implementation of the tested algorithm
void sequentialScan(int* array, int* out, int size)
{
	out[0] = 0;
	for(int i=1; i<size; i++)
	{
		out[i] = array[i-1] + out[i-1];
	}
}

std::chrono::steady_clock::time_point begin;
std::chrono::steady_clock::time_point end;

int main()
{
	algebra::compute::Context::initContext();

	for(int i=0; i<N; i++)
	{
		std::cout << "Size: " << sizes[i] << std::endl;
		srand((unsigned)time(0));

		int* in = (int*)malloc(sizeof(int) * sizes[i]);

		for(int j=0; j<sizes[i]; j++){
        in[j] = (rand()%10)+1;
		}

		int* out0 = (int*)malloc(sizeof(int) * sizes[i]);
		int* out1 = (int*)malloc(sizeof(int) * sizes[i]);

		// Do the GPU parallel scan
		begin = std::chrono::steady_clock::now();
		scan(in, out0, sizes[i]);
		end = std::chrono::steady_clock::now();

		std::cout << "Total Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;

		// Do the sequential scan
		begin = std::chrono::steady_clock::now();
		sequentialScan(in, out1, sizes[i]);
		end = std::chrono::steady_clock::now();

		std::cout << "Naive Total: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;
		std::cout << std::endl;
	
		delete in;
		delete out0;
		delete out1;
	}

	algebra::compute::Context::stopContext();
	return 0;
}
