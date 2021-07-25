#include <assert.h>

#include "algebra/matrix/MatrixSubtractor.hpp"

using namespace karu;
using namespace algebra;

void MatrixSubtractor::sub(MatrixData* C, const MatrixData* const A, const MatrixData* const B, bool A_T, bool B_T)
{
	assert(A->columns() == B->columns());
	assert(A->lines() == B->lines());

	for(i32 i=0; i<C->m_lines; i+= C->blockHeight()) 
	{
		for(i32 j=0; j<C->m_columns; j+= C->blockWidth()) 
		{
			u32 row_margin = C->m_lines - i;
			u32 col_margin = C->m_columns - j;
	
			i32 width = std::min(C->m_block_width, col_margin);
			i32 heigth = std::min(C->m_block_heigth, row_margin);

			for(i32 y=0; y<heigth; y++)
			{
				for(i32 x=0; x<width; x++)
				{
					i32 Ci = i+y;
					i32 Cj = j+x;

					C->set(
						Ci, Cj,
						A->get((1 - A_T) * Ci + A_T * Cj, (1 - A_T) * Cj + A_T * Ci) -
						B->get((1 - B_T) * Ci + B_T * Cj, (1 - B_T) * Cj + B_T * Ci)
					);
				}
			}
		}
	}
}
