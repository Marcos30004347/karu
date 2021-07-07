#pragma once

#include "algebra/core/types.hpp"

namespace karu {
namespace algebra {

class MatrixData {
	public:

	MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth, f32* data);
	MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth);

	f32 get(i32 i, i32 j) const;

	const u32 lines() const;
	const u32 columns() const;
	const u32 blockHeight() const;
	const u32 blockWidth() const;

	const void set(i32 i, i32 j, f32 val) const;

	static u32 getBlockStartIndex(u32 i, u32 j, u32 columns, u32 block_height, u32 block_width);
	static u32 getIndex(u32 i, u32 j, u32 columns, u32 block_height, u32 block_width);
	
	private:
	f32* m_data;

	const u32 m_lines;
	const u32 m_columns;
	const u32 m_block_heigth;
	const u32 m_block_width;
};

}
}



