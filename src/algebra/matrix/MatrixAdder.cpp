#include <assert.h>

#include "algebra/matrix/MatrixAdder.hpp"

using namespace karu;
using namespace algebra;

MatrixData MatrixAdder::add(MatrixData* A, MatrixData* B, bool A_T, bool B_T)
{
	assert(A->columns() == B->columns());
	assert(A->lines() == B->lines());

	assert(A->blockHeight() == B->blockHeight());
	assert(A->blockWidth() == B->blockWidth());

	MatrixData C = MatrixData(A->lines(), B->columns(), A->blockWidth(), B->blockHeight());

	// Iterate over C blocks
	for(i32 i=0; i<A->lines()/A->blockHeight(); i++) {
		for(i32 j=0; j<B->columns()/B->blockWidth(); j++) {

			for(i32 x=0; x<B->blockWidth(); x++) {
				for(i32 y=0; y<A->blockHeight(); y++) {

					i32 Ci = (A->blockHeight()*i) % A->lines() + y;
					i32 Cj = (B->blockWidth()*j) % B->columns() + x;

					C.set(Ci, Cj, A->get((1 - A_T) * Ci + A_T * Cj, (1 - A_T) * Cj + A_T * Ci) + B->get((1 - B_T) * Ci + B_T * Ci, (1 - B_T) * Cj + B_T * Ci));
				}
			}
		}
	}
	return C;
}
