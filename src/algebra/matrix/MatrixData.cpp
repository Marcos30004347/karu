#include "MatrixData.hpp"
#include <algorithm>
#include <cmath>

using namespace karu;
using namespace algebra;

MatrixData::MatrixData()
{
	this->m_data = nullptr;
	this->m_lines = 0;
	this->m_columns = 0;
	this->m_block_heigth = 0;
	this->m_block_width = 0;
	this->m_stored_lines  = 0;
	this->m_stored_column  = 0;
}

MatrixData::MatrixData(const MatrixData& other)
{

	this->m_lines = other.m_lines;
	this->m_columns = other.m_columns;
	this->m_block_heigth = other.m_block_heigth;
	this->m_block_width = other.m_block_width;
	this->m_data = new f32[this->m_stored_lines*this->m_stored_column];
	this->m_stored_lines  = other.m_stored_lines;
	this->m_stored_column  = other.m_stored_column;
	
	std::copy(other.m_data, other.m_data + (other.m_stored_lines*other.m_stored_column), this->m_data);
}
MatrixData::MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth, std::initializer_list<f32> data)
{
	m_lines = lines;
	m_columns = columns;
	m_block_heigth = block_heigth;
	m_block_width = block_width;
	this->m_stored_lines  = std::ceil(this->m_lines/(f32)this->m_block_heigth) * this->m_block_heigth;
	this->m_stored_column = std::ceil(this->m_columns/(f32)this->m_block_width) * this->m_block_width;
	
	this->m_data = new f32[this->m_stored_column*this->m_stored_lines];
	std::fill(this->m_data, this->m_data + this->m_stored_column*this->m_stored_lines, 0);
	
	// define extra space as identity
	for(int i=0; i<std::min(this->m_stored_lines - this->m_lines, this->m_stored_column - this->m_columns); i++)
	{
		i32 block_y = (this->m_lines + i)/this->m_block_heigth;
		i32 block_x = (this->m_columns + i)/this->m_block_width;

		i32 y = (this->m_lines + i) - block_y*this->m_block_heigth;
		i32 x = (this->m_columns + i) - block_x*this->m_block_width;

		this->m_data[this->stride(block_y, block_x) + y*this->m_block_width + x] = 1.f;
	}

	i32 idx = 0;

	for(i32 block_y=0; block_y<this->m_stored_lines/this->m_block_heigth; block_y++)
	{
		for(i32 block_x=0; block_x<this->m_stored_column/this->m_block_width; block_x++)
		{
			// Current block width and current block heigth
			u32 row_margin = lines - block_y*this->m_block_heigth;
			u32 col_margin = columns - block_x*this->m_block_width;
	
			i32 width = std::min(m_block_width, col_margin);
			i32 heigth = std::min(m_block_heigth, row_margin);

			// where current block starts 

			u32 idx = this->m_block_heigth * this->m_columns * block_y + block_x*this->m_block_width;

			for(i32 y=0; y<heigth; y++)
			{
				for(i32 x=0; x<width; x++)
				{
					this->m_data[this->stride(block_y, block_x) + y*this->m_block_width + x] = *(data.begin() + idx++);
				}
				//this->m_stored_lines - this->m_lines + m_block_width - width;
				idx += this->m_columns - (block_x+1)*this->m_block_width + (block_x)*this->m_block_width + m_block_width - width;
			}
		}
	}
}

MatrixData::MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth, f32* data)
{
	m_lines = lines;
	m_columns = columns;
	m_block_heigth = block_heigth;
	m_block_width = block_width;

	this->m_stored_lines  = std::ceil(this->m_lines/(f32)this->m_block_heigth) * this->m_block_heigth;
	this->m_stored_column = std::ceil(this->m_columns/(f32)this->m_block_width) * this->m_block_width;
	
	this->m_data = new f32[this->m_stored_lines*this->m_stored_column];
	std::fill(this->m_data, this->m_data + this->m_stored_column*this->m_stored_lines, 0);
	
	// define extra space as identity
	for(int i=0; i<std::min(this->m_stored_lines - this->m_lines, this->m_stored_column - this->m_columns); i++)
	{
		i32 block_y = (this->m_lines + i)/this->m_block_heigth;
		i32 block_x = (this->m_columns + i)/this->m_block_width;

		i32 y = (this->m_lines + i) - block_y*this->m_block_heigth;
		i32 x = (this->m_columns + i) - block_x*this->m_block_width;

		this->m_data[this->stride(block_y, block_x) + y*this->m_block_width + x] = 1.f;
	}

	i32 idx = 0;

	for(i32 block_y=0; block_y<this->m_stored_lines/this->m_block_heigth; block_y++)
	{
		for(i32 block_x=0; block_x<this->m_stored_column/this->m_block_width; block_x++)
		{
			// Current block width and current block heigth
			u32 row_margin = lines - block_y*this->m_block_heigth;
			u32 col_margin = columns - block_x*this->m_block_width;
	
			i32 width = std::min(m_block_width, col_margin);
			i32 heigth = std::min(m_block_heigth, row_margin);

			// u32 idx = block_y*this->m_block_heigth*this->m_block_width + block_x*this->m_block_width;
			u32 idx = this->m_block_heigth * this->m_columns * block_y + block_x*this->m_block_width;

			for(i32 y=0; y<heigth; y++)
			{
				for(i32 x=0; x<width; x++)
				{
					this->m_data[this->stride(block_y, block_x) + y*this->m_block_width + x] = data[idx];
				}
				idx += this->m_columns - (block_x+1)*this->m_block_width + (block_x)*this->m_block_width + m_block_width - width;
			}
		}
	}
}

MatrixData::MatrixData(u32 lines, u32 columns, u32 block_width, u32 block_heigth)
{
	m_lines = lines;
	m_columns = columns;
	m_block_heigth = block_heigth;
	m_block_width = block_width;

	this->m_stored_lines  = std::ceil(this->m_lines/(f32)this->m_block_heigth) * this->m_block_heigth;
	this->m_stored_column = std::ceil(this->m_columns/(f32)this->m_block_width) * this->m_block_width;	
	
	this->m_data = new f32[this->m_stored_lines*this->m_stored_column];
	std::fill(this->m_data, this->m_data + this->m_stored_column*this->m_stored_lines, 0);

	// define extra space as identity
	for(int i=0; i<std::min(this->m_stored_lines - this->m_lines, this->m_stored_column - this->m_columns); i++)
	{
		i32 block_y = (this->m_lines + i)/this->m_block_heigth;
		i32 block_x = (this->m_columns + i)/this->m_block_width;

		i32 y = (this->m_lines + i) - block_y*this->m_block_heigth;
		i32 x = (this->m_columns + i) - block_x*this->m_block_width;

		this->m_data[this->stride(block_y, block_x) + y*this->m_block_width + x] = 1.f;
	}
}

MatrixData::~MatrixData()
{

	if(this->m_data)
		delete[] this->m_data;
}

const u32 MatrixData::lines() const
{
	return this->m_lines;
}

const u32 MatrixData::columns() const
{
	return this->m_columns;
}

const u32 MatrixData::blockWidth() const
{
	return this->m_block_width;
}

const u32 MatrixData::blockHeight() const
{
	return this->m_block_heigth;
}

#include <string.h>

MatrixData& MatrixData::operator=(const MatrixData& other)
{
	this->m_stored_column = other.m_stored_column;
	this->m_stored_lines  = other.m_stored_lines;
	this->m_lines 				= other.m_lines;
	this->m_columns 			= other.m_columns;
	this->m_block_width 	= other.m_block_width;
	this->m_block_heigth 	= other.m_block_heigth;

	if(this->m_data) delete[] this->m_data;
	
	this->m_data = new f32[this->m_stored_lines*this->m_stored_column];

	std::copy(other.m_data, other.m_data + (other.m_stored_lines * other.m_stored_column), this->m_data);

	return *this;
}
