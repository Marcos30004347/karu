#include <assert.h>

#include "algebra/lib/Scan.hpp"
#include "algebra/core/compute/Compute.hpp"

using namespace karu;

int in[] = {1, 2, 0, 3, 0, 1, 1, 0, 3, 3, 3, 2, 1, 2, 2, 0, 2, 0, 0, 2};
int out[20];

int main()
{
	algebra::compute::Context::initContext();

	scan(in, out, 20);

	algebra::compute::Context::stopContext();
	return 0;
}
