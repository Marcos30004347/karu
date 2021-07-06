#include "Matrix.hpp"

using namespace karu;

Matrix::Matrix(u32 lines, u32 columns, f32* data):
	m_lines(lines),
	m_columns(columns),
	m_block_heigth(4),
	m_block_width(4)
{
	this->m_data = new f32[lines*columns];
	i32 idx = 0;

	for(i32 block_line=0; block_line<lines; block_line+=m_block_heigth)
	{
		for(i32 block_column=0; block_column<columns; block_column+=m_block_width)
		{
			for(i32 y=0; y<m_block_heigth; y++)
			{
				for(i32 x=0; x<m_block_width; x++)
				{
					this->m_data[idx++] = data[(block_line + y) * columns + block_column + x];
				}
			}
		}
	}

}
