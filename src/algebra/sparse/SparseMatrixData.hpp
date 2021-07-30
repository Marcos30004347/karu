#pragma once

#include <vector>
#include "algebra/core/types.hpp"
#include "algebra/compute/Buffer.hpp"

namespace karu {
namespace algebra {



class SparseMatrixData {
	public:
	SparseMatrixData(
		u64 block_width, u64 block_heigth,
		u64 lines, u64 columns,
		
		std::vector<u64> row_ptr,
		std::vector<u64> col_idx,
		std::vector<f32> data
	);

	f32 get(u64 l, u64 c);

	void print();

	inline u64 lines()
	{
		return this->bcsr_lines;
	}

	inline u64 columns()
	{
		return this->bcsr_columns;
	} 

	inline f32* data()
	{
		return this->bcsr_data.data();
	}

	private:
	u64 bcsr_block_heigth;
	u64 bcsr_block_width;

	u64 bcsr_lines;
	u64 bcsr_columns;

	// row_ptr[i] stores where each row first block for bcsr_col_idx,
	// so the i'th row starts at bcsr_col_idx[row_ptr[i]] columns and
	// have bcsr_col_idx[row_ptr[i+1]] - bcsr_col_idx[row_ptr[i]] columns.
	std::vector<u64> bcsr_row_ptr;

	// bcsr_col_idx[i] holds the column id for the i'th block
	std::vector<u64> bcsr_col_idx;

	// bcsr_data is the array storing the data of the matrix
	// data is stored in column major
	std::vector<f32> bcsr_data;

	friend class SparseMatrixMultiplayer;
};

}
}
