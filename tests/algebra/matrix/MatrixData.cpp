#include <assert.h>
#include <iostream>

#include "algebra/matrix/MatrixData.hpp"

using namespace karu;

int main()
{

	f32 A_data[] = {
		1.f, 	2.f, 	3.f, 	4.f,
		5.f, 	6.f, 	7.f, 	8.f,
		9.f, 	10.f, 11.f, 12.f,
		13.f, 14.f, 15.f, 16.f,
	};

	algebra::MatrixData data = algebra::MatrixData(4, 4, 4, 4, A_data);
	
	assert(data.get(0, 0) == 1.f);
	assert(data.get(0, 1) == 2.f);
	assert(data.get(0, 2) == 3.f);
	assert(data.get(0, 3) == 4.f);
	assert(data.get(1, 0) == 5.f);
	assert(data.get(1, 1) == 6.f);
	assert(data.get(1, 2) == 7.f);
	assert(data.get(1, 3) == 8.f);
	assert(data.get(2, 0) == 9.f);
	assert(data.get(2, 1) == 10.f);
	assert(data.get(2, 2) == 11.f);
	assert(data.get(2, 3) == 12.f);
	assert(data.get(3, 0) == 13.f);
	assert(data.get(3, 1) == 14.f);
	assert(data.get(3, 2) == 15.f);
	assert(data.get(3, 3) == 16.f);

	return 0;
}
