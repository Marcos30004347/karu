#pragma once

#include "algebra/core/types.hpp"
#include <initializer_list>
#include <iostream>

namespace karu {
namespace algebra {

class MatrixData {
	friend class MatrixAdder;
	friend class MatrixSubtractor;
	friend class MatrixMultiplayer;

	public:
	MatrixData();
	MatrixData(const MatrixData& other);
	MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth, std::initializer_list<f32> data);
	MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth, f32* data);
	MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth);
	~MatrixData();

	const u32 lines() const;
	const u32 columns() const;

	const u32 blockHeight() const;
	const u32 blockWidth() const;

	inline f32& get(i32 i, i32 j) const
	{
		i32 block_y = i/this->m_block_heigth;
		i32 block_x = j/this->m_block_width;

		i32 y = i - block_y*this->m_block_heigth;
		i32 x = j - block_x*this->m_block_width;
	
		return this->m_data[this->stride(block_y, block_x) + y*this->m_block_width + x];
	}

	inline const void set(i32 i, i32 j, f32 val) const
	{
		i32 block_y = i/this->m_block_heigth;
		i32 block_x = j/this->m_block_width;

		i32 y = i - block_y*this->m_block_heigth;
		i32 x = j - block_x*this->m_block_width;

	 	this->m_data[this->stride(block_y, block_x) + y*this->m_block_width + x] = val;
	}
	
	inline u32 stride(u32 block_line, u32 block_column) const
	{
		u32 blocks_per_block_line = this->m_stored_column/this->m_block_width;
		u32 block_size = this->m_block_heigth * this->m_block_width;
		return block_line*blocks_per_block_line*block_size + block_column*block_size;
	}
	// inline static u32 getBlockStartIndex(u32 i, u32 j, u32 columns, u32 block_height, u32 block_width)
	// {
	// 	u32 block_line = i/block_height;
	// 	u32 block_column = j/block_width;
	// 	return block_line*(columns/block_width)*(block_height * block_width) + block_column * (block_height * block_width);
	// }

	// inline static u32 getIndex(u32 i, u32 j, u32 columns, u32 block_height, u32 block_width)
	// {
	// 	u32 y = i%block_height;
	// 	u32 x = j%block_width;
	// 	return getBlockStartIndex(i,j, columns, block_height, block_width) + (y*block_width + x);
	// }
	
	MatrixData& operator=(const MatrixData& other);

	// private:
	f32* m_data;

	u32 m_lines;
	u32 m_columns;
	u32 m_stored_lines;
	u32 m_stored_column;
	u32 m_block_heigth;
	u32 m_block_width;
};

}
}



