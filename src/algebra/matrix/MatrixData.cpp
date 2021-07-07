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

u32 MatrixData::getBlockStartIndex(u32 i, u32 j, u32 columns, u32 block_height, u32 block_width) {
	u32 block_line = i/block_height;
	u32 block_column = j/block_width;
	return block_line*(columns/block_width)*(block_height * block_width) + block_column * (block_height * block_width);
}

u32 MatrixData::getIndex(u32 i, u32 j, u32 columns, u32 block_height, u32 block_width) {
	u32 y = i%block_height;
	u32 x = j%block_width;
	return getBlockStartIndex(i,j, columns, block_height, block_width) + (y*block_width + x);
}

const void MatrixData::set(i32 i, i32 j, f32 val) const
{
	this->m_data[MatrixData::getIndex(i,j, this->columns(), this->blockHeight(),  this->blockWidth())] = val;
}

f32 MatrixData::get(i32 i, i32 j) const
{
	return this->m_data[this->getIndex(i, j, this->columns(), this->blockHeight(), this->blockWidth())];
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
