#include "algebra/matrix/MatrixTransposer.hpp"

namespace karu::algebra {

void MatrixTransposer::transpose(MatrixData* C, const MatrixData* const A)
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
					C->set(i+y, j+x, A->get(j+x, i+y));
				}
			}
		}
	}
}

}
