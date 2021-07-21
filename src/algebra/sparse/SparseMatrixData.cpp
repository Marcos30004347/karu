#include <iostream>

#include "algebra/sparse/SparseMatrixData.hpp"

using namespace karu;
using namespace algebra;

SparseMatrixData::SparseMatrixData(
	u64 block_width, u64 block_heigth,
	u64 lines, u64 columns,
	std::vector<u64> row_ptr,
	std::vector<u64> col_idx,
	std::vector<f32> data
)
{
	this->bcsr_block_width = block_width;
	this->bcsr_block_heigth = block_heigth;
	this->bcsr_row_ptr = row_ptr;
	this->bcsr_col_idx = col_idx;
	this->bcsr_data = data;

	this->bcsr_lines = lines;
	this->bcsr_columns = columns;


}

f32 SparseMatrixData::get(u64 l, u64 c)
{
	u64 i = l/this->bcsr_block_heigth;

	u64 first_block = this->bcsr_row_ptr[i]; 
	u64 last_block = this->bcsr_row_ptr[i+1];

	u64 n = last_block - first_block;

	bool block_column_found = false;

	u64 t = 0;
	u64 j = c;
	
	while(
		j > this->bcsr_col_idx[first_block + t] + this->bcsr_block_heigth - 1 &&
		t < last_block - first_block
	) t++;
	
	if(j < this->bcsr_col_idx[first_block + t] || j > this->bcsr_col_idx[first_block + t] + this->bcsr_block_width - 1)
		return 0;

	u64 b = first_block + t;

	u64 bl = l%(this->bcsr_block_heigth);
	u64 bc = (c - this->bcsr_col_idx[first_block + t])%(this->bcsr_block_width);

	u64 bs = this->bcsr_block_heigth * this->bcsr_block_width;

	return this->bcsr_data[b*bs + bc*this->bcsr_block_heigth + bl];
}

void SparseMatrixData::print()
{
	for(u32 l = 0; l < this->lines()/this->bcsr_block_heigth; l++)
	{
		for(u32 bl = 0; bl < this->bcsr_block_heigth; bl++)
		{
			if(this->bcsr_row_ptr[l+1] - this->bcsr_row_ptr[l] == 0)
			{
				for(u32 c = 0; c < this->columns(); c++)
				{
					std::cout << "0 ";
				}
			}
			else
			{
				u32 col = 0;
				for(u32 b = this->bcsr_row_ptr[l]; b < this->bcsr_row_ptr[l+1]; b++)
				{

					while(col < this->bcsr_col_idx[b])
					{
						std::cout << "0 ";
						col++;
					}

					for(i32 c = 0; c<this->bcsr_block_width; c++)
					{
						std::cout << this->bcsr_data[
							b*this->bcsr_block_heigth*bcsr_block_width +
							c*bcsr_block_heigth + bl
						] << " ";
						col++;
					}
				
				}
				while(col < this->columns())
				{
					std::cout << "0 ";
					col++;
				}
			}
			std::cout << std::endl;
		}
	}
}
