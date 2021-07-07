#include <assert.h>

#include "algebra/matrix/MatrixMultiplayer.hpp"

using namespace karu;
using namespace algebra;

MatrixData MatrixMultiplayer::multiply(MatrixData* A, MatrixData* B, bool A_T, bool B_T)
{
	assert(A->columns() == B->lines());
	assert(A->blockWidth() == B->blockHeight());
	
	MatrixData C = MatrixData(A->lines(), B->columns(), A->blockWidth(), B->blockHeight());

	// Iterate over C blocks
	for(i32 i=0; i<A->lines()/A->blockHeight(); i++) {
		for(i32 j=0; j<B->columns()/B->blockWidth(); j++) {
			//working C(i,j) block
			// C(i,j) = A(i,:)*B(:,j)
			for(i32 k=0; k<A->columns()/A->blockWidth(); k++) {
				// k'th A block line
				// k'th B block column

				// Multiply block matrix
				for(i32 y=0; y<B->blockHeight(); y++) {
					for(i32 x=0; x<A->blockWidth(); x++) {
						for(i32 q=0; q<A->blockWidth(); q++) {
							i32 Ai = (B->blockHeight()*i) % A->lines() + y;
							i32 Aj = (A->blockWidth()*k) % A->columns() + (x+q) % A->blockWidth();

							i32 Bi = (B->blockHeight()*k) % B->lines() + (x+q) % A->blockWidth();
							i32 Bj = (A->blockWidth()*j) % B->columns() + x;
							
							i32 Ci = (B->blockHeight()*i) % A->lines() + y;
							i32 Cj = (A->blockWidth()*j) % B->columns() + x;

							C.set(
								Ci, Cj,
								C.get(Ci,Cj) +
								A->get((1 - A_T) * Ai + A_T * Aj, (1 - A_T) * Aj + A_T * Ai) *
								B->get((1 - B_T) * Bi + B_T * Bj, (1 - B_T) * Bj + B_T * Bi)
							);
						}
					}
				}
			}
		}
	}	
	return C;
}
