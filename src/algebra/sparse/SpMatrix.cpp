#include "SpMatrix.hpp"
#include "SparseMatrixMultiplayer.hpp"

namespace karu::algebra {

SpMatrix::SpMatrix(){}

SpMatrix::SpMatrix(
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


Matrix SpMatrix::operator*(Matrix& other)
{
	if(other.columns() == 1)
	{
		Matrix y = Matrix(this->m_data.lines(), other.columns(), other.m_data.blockWidth(), other.m_data.blockHeight());

		SparseMatrixMultiplayer::sparseMVMultiplyCPU(
			&this->m_data,
			&other.m_data,
			&y.m_data
		);
		return y;
	}
	else
	{
		std::cout << "ERROR: SpMatrix * Matrix not implemented\n";
		exit(EXIT_FAILURE);
	}
}

void printMatrix(SpMatrix& mat)
{
	mat.m_data.print();
}


}
