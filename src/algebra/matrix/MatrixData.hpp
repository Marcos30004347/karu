#pragma once

#include "algebra/core/types.hpp"

namespace karu {
namespace algebra {

class MatrixData {
	public:

	MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth, f32* data);
	MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth);
	~MatrixData();

	const u32 lines() const;
	const u32 columns() const;
	const u32 blockHeight() const;
	const u32 blockWidth() const;

	inline f32 get(i32 i, i32 j) const
	{
		return this->m_data[this->getIndex(i, j, this->columns(), this->blockHeight(), this->blockWidth())];
	}

	inline const void set(i32 i, i32 j, f32 val) const
	{
		this->m_data[MatrixData::getIndex(i,j, this->columns(), this->blockHeight(),  this->blockWidth())] = val;
	}

	inline static u32 getBlockStartIndex(u32 i, u32 j, u32 columns, u32 block_height, u32 block_width)
	{
		u32 block_line = i/block_height;
		u32 block_column = j/block_width;
		return block_line*(columns/block_width)*(block_height * block_width) + block_column * (block_height * block_width);
	}

	inline static u32 getIndex(u32 i, u32 j, u32 columns, u32 block_height, u32 block_width)
	{
		u32 y = i%block_height;
		u32 x = j%block_width;
		return getBlockStartIndex(i,j, columns, block_height, block_width) + (y*block_width + x);
	}
	
	private:

	f32* m_data;
	const u32 m_lines;
	const u32 m_columns;
	const u32 m_block_heigth;
	const u32 m_block_width;
};

}
}



