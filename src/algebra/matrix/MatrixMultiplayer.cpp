#include <assert.h>

#include "algebra/matrix/MatrixMultiplayer.hpp"

using namespace karu;
using namespace algebra;

void MatrixMultiplayer::mul(MatrixData* C, const MatrixData* const A, const MatrixData* const B, bool A_T, bool B_T)
{
	assert(A->columns() == B->lines());

	// code for GPU multiplication
	for(i32 i=0; i<C->m_lines; i+= C->blockHeight()) 
	{
		for(i32 j=0; j<C->m_columns; j+= C->blockWidth()) 
		{
			// loops above will become one block
			for(i32 y=0; y<C->m_block_heigth; y++)
			{
				for(i32 x=0; x<C->m_block_width; x++)
				{
					// loops above will become one thread
					f32 acc = 0;

					i32 Ci = i + y, Cj = j + x;
					for(i32 k=0; k < A->columns(); k += std::min(A->blockWidth(), B->blockHeight()) /*C->blockWidth()*/)
					{
						for(i32 q=0; q < std::min(A->blockWidth(), B->blockHeight()); q++)
						{
							i32 Ai = (i + y), Aj = (k + q);
							i32 Bi = (k + q), Bj = (j + x);

							bool Ai_bounds = Ai < A->m_stored_lines;
							bool Aj_bounds = Aj < A->m_stored_column;
							bool Bi_bounds = Bi < B->m_stored_lines;
							bool Bj_bounds = Bj < B->m_stored_column;

							Ai = (A->m_stored_lines-1)*!Ai_bounds  + Ai*Ai_bounds;
							Aj = (A->m_stored_column-1)*!Aj_bounds + Aj*Aj_bounds;
							Bi = (B->m_stored_lines-1)*!Bi_bounds  + Bi*Bi_bounds;
							Bj = (B->m_stored_column-1)*!Bj_bounds + Bj*Bj_bounds;

							acc += A->get(Ai, Aj) * B->get(Bi, Bj);
							// i32 Ai = (i + y)%A->m_stored_lines, Aj = (k + q)%A->m_stored_column;
							// i32 Bi = (k + q)%B->m_stored_lines, Bj = (j + x)%B->m_stored_column;
							
							// acc += A->get(Ai, Aj) * B->get(Bi, Bj);
						}
					}
					C->set(Ci, Cj, acc);
				}
			}
		}
	}


	// Above lines just set the padding of the matrix to an identity like matrix
	for(int j=0; j<C->m_stored_column - C->m_columns; j++)
	{
		for(int i=0; i<C->m_lines; i++)
		{
			i32 block_y = i/C->m_block_heigth;
			i32 block_x = (C->m_columns + j)/C->m_block_width;

			i32 y = i - block_y*C->m_block_heigth;
			i32 x = (C->m_columns + j) - block_x*C->m_block_width;

			C->m_data[C->stride(block_y, block_x) + y*C->m_block_width + x] = 0.f;
		}
	}
	for(int i=0; i<C->m_stored_lines - C->m_lines; i++)
	{
		for(int j=0; j<C->m_columns; j++)
		{
			i32 block_y = (C->m_lines + i)/C->m_block_heigth;
			i32 block_x = j/C->m_block_width;

			i32 y = (C->m_lines + i) - block_y*C->m_block_heigth;
			i32 x = j - block_x*C->m_block_width;

			C->m_data[C->stride(block_y, block_x) + y*C->m_block_width + x] = 0.f;
		}
	}
	for(int i=0; i<std::min(C->m_stored_lines - C->m_lines, C->m_stored_column - C->m_columns); i++)
	{
		i32 block_y = (C->m_lines + i)/C->m_block_heigth;
		i32 block_x = (C->m_columns + i)/C->m_block_width;

		i32 y = (C->m_lines + i) - block_y*C->m_block_heigth;
		i32 x = (C->m_columns + i) - block_x*C->m_block_width;

		C->m_data[C->stride(block_y, block_x) + y*C->m_block_width + x] = 1.f;
	}

	return;

	// // code for CPU multiplication
	// for(i32 i=0; i<C->m_lines; i+= C->blockHeight()) 
	// {
	// 	for(i32 j=0; j<C->m_columns; j+= C->blockWidth()) 
	// 	{
	// 		for(i32 k=0; k<A->columns(); k+=A->blockWidth())
	// 		{
	// 			u32 row_margin = C->m_lines - i;
	// 			u32 col_margin = C->m_columns - j;
		
	// 			i32 width = std::min(A->m_block_width, col_margin);
	// 			i32 heigth = std::min(B->m_block_heigth, row_margin);

	// 			for(i32 y=0; y<heigth; y++)
	// 			{
	// 				for(i32 x=0; x<width; x++)
	// 				{
	// 					for(i32 q=0; q<A->m_block_width; q++)
	// 					{
	// 						i32 Ci = i + y, Cj = j + x;
	// 						i32 Ai = i + y, Aj = k + q;
	// 						i32 Bi = k + q, Bj = j + x;

	// 						C->set(
	// 							Ci, Cj,
	// 							C->get(Ci,Cj) +
	// 							A->get((1 - A_T) * Ai + A_T * Aj, (1 - A_T) * Aj + A_T * Ai) *
	// 							B->get((1 - B_T) * Bi + B_T * Bj, (1 - B_T) * Bj + B_T * Bi)
	// 						);
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	// for(int i=0; i<std::min(C->m_stored_lines - C->m_lines, C->m_stored_column - C->m_columns); i++)
	// {
	// 	i32 block_y = (C->m_lines + i)/C->m_block_heigth;
	// 	i32 block_x = (C->m_columns + i)/C->m_block_width;

	// 	i32 y = (C->m_lines + i) - block_y*C->m_block_heigth;
	// 	i32 x = (C->m_columns + i) - block_x*C->m_block_width;

	// 	C->m_data[C->stride(block_y, block_x) + y*C->m_block_width + x] = 1.f;
	// }

}

void MatrixMultiplayer::mul(MatrixData* C, const MatrixData* const A, const f32 alpha)
{
	// code for GPU multiplication
	for(i32 i=0; i<C->m_lines; i+= C->blockHeight()) 
	{
		for(i32 j=0; j<C->m_columns; j+= C->blockWidth()) 
		{
			// loops above will become one block
			for(i32 y=0; y<C->m_block_heigth; y++)
			{
				for(i32 x=0; x<C->m_block_width; x++)
				{
					// loops above will become one thread
					i32 Ci = i + y, Cj = j + x;
					C->set(Ci, Cj, A->get(Ci, Cj)*alpha);
				}
			}
		}
	}

	for(int i=0; i<std::min(C->m_stored_lines - C->m_lines, C->m_stored_column - C->m_columns); i++)
	{
		i32 block_y = (C->m_lines + i)/C->m_block_heigth;
		i32 block_x = (C->m_columns + i)/C->m_block_width;

		i32 y = (C->m_lines + i) - block_y*C->m_block_heigth;
		i32 x = (C->m_columns + i) - block_x*C->m_block_width;

		C->m_data[C->stride(block_y, block_x) + y*C->m_block_width + x] = 1.f;
	}

	return;
}
