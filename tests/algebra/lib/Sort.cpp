#include <assert.h>

#include "algebra/lib/Sort.hpp"
#include "algebra/core/compute/Compute.hpp"

using namespace karu;

int keys[] = {1, 2, 0, 3, 0, 1, 1, 0, 3, 3, 3, 2, 1, 2, 2, 0, 2, 0, 0, 2};
int vals[] = {1, 2, 0, 3, 0, 1, 1, 0, 3, 3, 3, 2, 1, 2, 2, 0, 2, 0, 0, 2};

int main()
{
	algebra::compute::Context::initContext();

	sort(keys, vals, 20, 32);

	algebra::compute::Context::stopContext();
	return 0;
}
