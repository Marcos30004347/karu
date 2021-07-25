#include <assert.h>
#include <iostream>

#include "algebra/matrix/MatrixData.hpp"
#include "algebra/matrix/MatrixMultiplayer.hpp"

using namespace karu;

int main()
{
	// f32 data[] = {
	// 	1.f, 	2.f, 	3.f, 	4.f,
	// 	5.f, 	6.f, 	7.f, 	8.f,
	// 	9.f, 	10.f, 11.f, 12.f,
	// 	13.f, 14.f, 15.f, 16.f,
	// };

	// algebra::MatrixData A = algebra::MatrixData(4, 4, 4, 4, data);
	// algebra::MatrixData B = algebra::MatrixData(4, 4, 4, 4, data);
	
	// algebra::MatrixData C = algebra::MatrixMultiplayer::mul(&A, &B, false, false);

	// assert(C.get(0, 0) == 90);
	// assert(C.get(0, 1) == 100);
	// assert(C.get(0, 2) == 110);
	// assert(C.get(0, 3) == 120);

	// assert(C.get(1, 0) == 202);
	// assert(C.get(1, 1) == 228);
	// assert(C.get(1, 2) == 254);
	// assert(C.get(1, 3) == 280);

	// assert(C.get(2, 0) == 314);
	// assert(C.get(2, 1) == 356);
	// assert(C.get(2, 2) == 398);
	// assert(C.get(2, 3) == 440);
	
	// assert(C.get(3, 0) == 426);
	// assert(C.get(3, 1) == 484);
	// assert(C.get(3, 2) == 542);
	// assert(C.get(3, 3) == 600);

	// algebra::MatrixData D = algebra::MatrixMultiplayer::mul(&A, &B, true, true);

	// assert(D.get(0, 0) == 90);
	// assert(D.get(0, 1) == 202);
	// assert(D.get(0, 2) == 314);
	// assert(D.get(0, 3) == 426);

	// assert(D.get(1, 0) == 100);
	// assert(D.get(1, 1) == 228);
	// assert(D.get(1, 2) == 356);
	// assert(D.get(1, 3) == 484);

	// assert(D.get(2, 0) == 110);
	// assert(D.get(2, 1) == 254);
	// assert(D.get(2, 2) == 398);
	// assert(D.get(2, 3) == 542);
	
	// assert(D.get(3, 0) == 120);
	// assert(D.get(3, 1) == 280);
	// assert(D.get(3, 2) == 440);
	// assert(D.get(3, 3) == 600);

	// algebra::MatrixData E = algebra::MatrixMultiplayer::mul(&A, &B, false, true);

	// assert(E.get(0, 0) == 30);
	// assert(E.get(0, 1) == 70);
	// assert(E.get(0, 2) == 110);
	// assert(E.get(0, 3) == 150);

	// assert(E.get(1, 0) == 70);
	// assert(E.get(1, 1) == 174);
	// assert(E.get(1, 2) == 278);
	// assert(E.get(1, 3) == 382);

	// assert(E.get(2, 0) == 110);
	// assert(E.get(2, 1) == 278);
	// assert(E.get(2, 2) == 446);
	// assert(E.get(2, 3) == 614);
	
	// assert(E.get(3, 0) == 150);
	// assert(E.get(3, 1) == 382);
	// assert(E.get(3, 2) == 614);
	// assert(E.get(3, 3) == 846);

	return 0;
}
