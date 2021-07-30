#pragma once

#include "algebra/matrix/Matrix.hpp"
#include "algebra/sparse/SparseMatrixData.hpp"
#include "algebra/vector/Vector.hpp"
#include <initializer_list>

namespace karu::algebra {

class SMatrix {
public:
	SparseMatrixData m_data;

	SMatrix(
		u64 lines, u64 columns,
		u64 block_heigth, u64 block_width,
		std::vector<u64> row_ptr,
		std::vector<u64> col_idx,
		std::vector<f32> data
	);

	Matrix& operator*(const Matrix& other);
};

}
