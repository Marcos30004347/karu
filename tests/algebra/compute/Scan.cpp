#include <chrono>
#include <random>
#include <assert.h>
#include <iostream>

#include "algebra/compute/lib/Scan.hpp"
#include "algebra/compute/Compute.hpp"

using namespace karu::algebra::compute;

// Number of tests to execute
#define N 8

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
	Context::initContext();

	for(int i=0; i<N; i++)
	{
		std::cout << "Size: " << sizes[i] << std::endl;
		srand((unsigned)time(0));

		int* in = (int*)malloc(sizeof(int) * sizes[i]);

		for(int j=0; j<sizes[i]; j++){
        in[j] = (rand()%10)+1;
		}

		// Do the GPU parallel scan
		Buffer inBuff = Buffer(in, sizeof(int)*sizes[i], Buffer::READ_WRITE, false);
		Buffer outBuff = Buffer(sizeof(int)*sizes[i], Buffer::READ_WRITE, Buffer::MEM_GPU);
		inBuff.upload();
		begin = std::chrono::steady_clock::now();
		
		// Do the scan
		scan(&inBuff, &outBuff, sizes[i]);
	
		end = std::chrono::steady_clock::now();

		std::cout << "GPU Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;

		// Do the sequential scan
		int* cpuOut = (int*)malloc(sizeof(int) * sizes[i]);
		begin = std::chrono::steady_clock::now();
		sequentialScan(in, cpuOut, sizes[i]);
		end = std::chrono::steady_clock::now();


		std::cout << "CPU Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;
		std::cout << std::endl;
		
		int* gpuOut = (int*)outBuff.download();
		for(int j=0; j<sizes[i]; j++)
		{
			assert(gpuOut[j] == cpuOut[j]);
		}

		delete in;
		delete cpuOut;
		delete gpuOut;
	}

	Context::stopContext();
	return 0;
}
