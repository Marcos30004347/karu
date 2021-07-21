#include <assert.h>
#include <chrono>

#include "algebra/compute/lib/Reduce.hpp"
#include "algebra/compute/Compute.hpp"

using namespace karu;
using namespace karu::algebra::compute;

#define N 2

unsigned sizes[N] = {
	256,
	512,
};

unsigned cpuReduce(int* arr, int size)
{
	unsigned val = 0;
	for(int i=0; i<size; i++)
	{
		val += arr[i];
	}
	return val;
}

std::chrono::steady_clock::time_point begin;
std::chrono::steady_clock::time_point end;

int main()
{
	algebra::compute::Context::initContext();

	// reduce(keys, 20);
	for(int i=0; i<N; i++)
	{
		std::cout << "Size: " << sizes[i] << std::endl;
		srand((unsigned)time(0));

		int* in = (int*)malloc(sizeof(int) * sizes[i]);

		for(int j=0; j<sizes[i]; j++){
        in[j] = (rand()%10)+1;
		}
		
		Buffer inp = Buffer(in, sizeof(int)*sizes[i], Buffer::READ_ONLY, false);
		inp.upload();
		
		Buffer out = Buffer(sizeof(int), Buffer::READ_ONLY, Buffer::MEM_GPU);
	
		begin = std::chrono::steady_clock::now();
		reduce(&inp, &out, sizes[i]);
		end = std::chrono::steady_clock::now();
		std::cout << "Total Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;

		unsigned int gpu_reduce = *(int*)out.download();
		
		begin = std::chrono::steady_clock::now();
		unsigned int cpu_reduce = cpuReduce(in, sizes[i]);
		end = std::chrono::steady_clock::now();
		std::cout << "Naive Total: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;
		std::cout << std::endl;
	
		assert(gpu_reduce == cpu_reduce);
	
		delete in;
	}
	algebra::compute::Context::stopContext();
	return 0;
}
