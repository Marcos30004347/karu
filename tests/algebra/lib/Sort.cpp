#include <assert.h>

#include "algebra/lib/Sort.hpp"
#include "algebra/core/compute/Compute.hpp"

using namespace karu;

int keys_data[] = {1, 2, 0, 3, 0, 1, 1, 0, 3, 3, 3, 2, 1, 2, 2, 0, 2, 0, 0, 2};
int vals_data[] = {1, 2, 0, 3, 0, 1, 1, 0, 3, 3, 3, 2, 1, 2, 2, 0, 2, 0, 0, 2};

int order_checking[] = {1, 2, 0, 3, 0, 1, 1, 0, 3, 3, 3, 2, 1, 2, 2, 0, 2, 0, 0, 2};
int main()
{

	algebra::compute::Context::initContext();

	Buffer keys_buffer = Buffer(keys_data, sizeof(int)*20, algebra::compute::Buffer::READ_WRITE, false);
	Buffer vals_buffer = Buffer(vals_data, sizeof(int)*20, algebra::compute::Buffer::READ_WRITE, false);
	
	keys_buffer.upload();
	vals_buffer.upload();

	sort(&keys_buffer, &vals_buffer, 20);

	int* keys = (int*)keys_buffer.download();
	int* vals = (int*)vals_buffer.download();

	for(int i=0; i<20; i++)
	{
		std::cout << "d[" << keys[i] << "] = " << vals[i] << "\n";
	}

	algebra::compute::Context::stopContext();
	return 0;
}
