#pragma once

#include <vector>
#include "algebra/core/types.hpp"

namespace karu {
namespace algebra {

class BlockSparseMatrixData {
	public:
	BlockSparseMatrixData(
		u64 lines, u64 columns,
		std::vector<u64> row_ptr,
		std::vector<u64> col_idx,
		std::vector<f32> data
	);

	f32 get(u64 l, u64 c);
	
	private:
	u64 bcsr_block_heigth = 2;
	u64 bcsr_block_width = 3;

	u64 bcsr_lines;
	u64 bcsr_columns;

	// row_ptr[i] stores how many blocks are stored before row i
	// it can be used as index for bcsr_col_idx as bcsr_col_idx[row_ptr[i]]
	std::vector<u64> bcsr_row_ptr;

	// bcsr_col_idx[i] holds the column id for the i'th block
	std::vector<u64> bcsr_col_idx;

	// bcsr_data is the array storing the data of the matrix
	std::vector<f32> bcsr_data;

};

}
}
