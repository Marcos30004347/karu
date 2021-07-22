#include <random>
#include <chrono>
#include <assert.h>
#include <iostream>

#include "algebra/compute/lib/Sort.hpp"

using namespace karu::algebra::compute;

#define N 5

unsigned sizes[N] = {
	256,
	512,
	1024,
	10000,
	100000,
};

std::chrono::steady_clock::time_point begin;
std::chrono::steady_clock::time_point end;

int partition(int* keys, int* vals, int l, int r)
{
	int pivo = keys[r];

	int i = l-1;

	for(int j=l; j<=r-1;j++)
	{
		if(keys[j] < pivo)
		{
			i++;

			int t = keys[i];
			keys[i] = keys[j];
			keys[j] = t;

			t = vals[i];
			vals[i] = vals[j];
			vals[j] = t;
		}
	}

	int t = keys[i+1];
	keys[i+1] = keys[r];
	keys[r] = t;

	t = vals[i+1];
	vals[i+1] = vals[r];
	vals[r] = t;

	return (i+1);
}
void sort_rec(int* keys, int* vals, int l, int r)
{
	if(l < r)
	{
		int p = partition(keys, vals, l, r);
		sort_rec(keys, vals, l, p-1);
		sort_rec(keys, vals, p+1, r);
	}
}
void sort(int* keys, int* vals, int size)
{
	sort_rec(keys, vals, 0, size-1);
}

int main()
{
	Context::initContext();

	for(int i=0; i<N; i++)
	{
		std::cout <<"Size: " << sizes[i] <<std::endl;

		int* keys_GPU = (int*)malloc(sizeof(int) * sizes[i]);
		int* keys_CPU = (int*)malloc(sizeof(int) * sizes[i]);
		int* vals_GPU = (int*)malloc(sizeof(int) * sizes[i]);
		int* vals_CPU = (int*)malloc(sizeof(int) * sizes[i]);
		
		for(int j=0; j<sizes[i]; j++)
		{
        keys_GPU[j] = sizes[i] - j;
        keys_CPU[j] = sizes[i] - j;
        vals_GPU[j] = sizes[i] - j;
        vals_CPU[j] = sizes[i] - j;
		}
	
		Buffer keys_buffer = Buffer(keys_GPU, sizeof(int)*sizes[i], Buffer::READ_WRITE, false);
		Buffer vals_buffer = Buffer(vals_GPU, sizeof(int)*sizes[i], Buffer::READ_WRITE, false);
		
		keys_buffer.upload();
		vals_buffer.upload();
		
		begin = std::chrono::steady_clock::now();
		sort(&keys_buffer, &vals_buffer, sizes[i]);
		end = std::chrono::steady_clock::now();
	
		std::cout << "GPU Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;

		begin = std::chrono::steady_clock::now();
		sort(keys_CPU, vals_CPU, sizes[i]);
		end = std::chrono::steady_clock::now();
	
		std::cout << "CPU Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;

		int* keys = (int*)keys_buffer.download();
		int* vals = (int*)vals_buffer.download();

		for(int j=1; j<sizes[i]; j++)
		{
			assert(keys[j-1] <= keys[j]);
		}
		std::cout << std::endl;
		delete keys_GPU;
		delete keys_CPU;
		delete vals_GPU;
		delete vals_CPU;
	}


	Context::stopContext();
	return 0;
}
