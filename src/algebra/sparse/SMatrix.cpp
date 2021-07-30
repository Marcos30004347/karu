#include "SMatrix.hpp"
#include "SparseMatrixMultiplayer.hpp"

namespace karu::algebra {

SMatrix::SMatrix(
	u64 lines,
	u64 columns,
	u64 block_heigth,
	u64 block_width,
	std::vector<u64> row_ptr,
	std::vector<u64> col_idx,
	std::vector<f32> data
	) : m_data{
		block_width,
		block_heigth,
		lines,
		columns,
		row_ptr,
		col_idx,
		data
	} {}


Matrix& SMatrix::operator*(const Matrix& other)
{
	if(other.columns() == 1)
	{
		Matrix y = Matrix(other.rows(), other.columns(), other.m_data.blockWidth(), other.m_data.blockHeight());
		SparseMatrixMultiplayer::sparseMVMultiplyGPU(
			&this->m_data,
			(MatrixData*)&other.m_data,
			&y.m_data
		);
	}
	else
	{
		std::cout << "ERROR: SMatrix * Matrix not implemented\n";
		exit(EXIT_FAILURE);
	}
}

}
