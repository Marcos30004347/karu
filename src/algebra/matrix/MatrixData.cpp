#include "MatrixData.hpp"

using namespace karu;
using namespace algebra;

MatrixData::MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth, f32* data):
	m_lines(lines),
	m_columns(columns),
	m_block_heigth(block_heigth),
	m_block_width(block_width)
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

MatrixData::MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth):
	m_lines(lines),
	m_columns(columns),
	m_block_heigth(block_heigth),
	m_block_width(block_width)
{
	this->m_data = new f32[lines*columns];
}

MatrixData::~MatrixData()
{
	delete this->m_data;
}

const u32 MatrixData::lines() const
{
	return this->m_lines;
}

const u32 MatrixData::columns() const
{
	return this->m_columns;
}

const u32 MatrixData::blockWidth() const
{
	return this->m_block_width;
}

const u32 MatrixData::blockHeight() const
{
	return this->m_block_heigth;
}
