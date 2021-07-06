#pragma once

#include "core/types.hpp"

namespace karu {

class Matrix {
	public:

	Matrix(u32 lines, u32 columns, f32* data);

	f32 get(i32 i, i32 j) const;

	const i32 lines() const;
	const i32 columns() const;

	private:
	f32* m_data;

	const u32 m_lines;
	const u32 m_columns;
	const u32 m_block_heigth;
	const u32 m_block_width;
};

}



