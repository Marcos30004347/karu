#include <iostream>
#include <cstdio>
#include <string.h>
#include <assert.h>

#include "matrix.hpp"
#include "types.hpp"

using namespace karu;

const char* mat_multiply_kernel = R"(
int get_block_start_index(int i, int j, int columns, int block_height, int block_width) {
	int block_line = i/block_height;
	int block_column = j/block_width;
	return block_line*(columns/block_width)*(block_height * block_width) + block_column * (block_height * block_width);
}

int get_index(int i, int j, int columns, int block_height, int block_width) {
	int y = i%block_height;
	int x = j%block_width;
	return get_block_start_index(i,j, columns, block_height, block_width) + (y*block_width + x);
}

__kernel void matrix_mul(
	const int block_height,
	const int block_width,
	const int lines,
	const int columns,
	const int A_columns,
	const int B_lines,
  	const __global float* A,
	const __global float* B,
  	__global float* C
) {
	const int i = get_global_id(0);
	const int j = get_global_id(1);

	for(int k=0; k<A_columns/block_width; k++) {
	
		for(int y=0; y<block_height; y++) {
			for(int x=0; x<block_width; x++) {
				for(int q=0; q<block_width; q++) {
					C[
						get_index((block_height*i)%lines + y, (block_width*j)%columns + x, columns, block_height, block_width)
					] += A[
						get_index((block_height*i)%lines + y, (block_width*k)%A_columns + (x+q)%block_width, A_columns, block_height, block_width)
					] * B[
						get_index((block_height*k)%B_lines + (x+q)%block_width, (block_width*j)%columns + x, columns, block_height, block_width)
					];
				}
			}
		}
	}
}
)";

i32 i32max(i32 a, i32 b) {
  i32 diff = a - b;
  i32 dsgn = diff >> 31;
  return a - (diff & dsgn);
}

i32 i32min(i32 a, i32 b) {
  i32 diff = a - b;
  i32 dsgn = diff >> 31;
  return b + (diff & dsgn);
}

matrix matrix::zeros(const i32 l, const i32 c, const i32 b_w, const i32 b_h) {
	matrix C(l, c, b_w, b_h);
	C._data = new f32[l*c];

	for(i32 i=0; i<i32min(l, c);i++)
		C._data[C.index_of(i,i)] = 0.f;

	return C;
}

matrix::matrix(
	const i32 l,
	const i32 c,
	const i32 b_w,
	const i32 b_h
):	_lines(l),
	_columns(c),
	_block_height(b_h),
	_block_width(b_w),
	_is_transposed(0)
{
	this->_data = new f32[l*c];

	for(i32 i=0; i < l*c; i++)
		this->_data[i] = 0.f;

	for(i32 i=0; i < i32min(l,c); i++)
		this->_data[index_of(i,i)] = 1.f;
}


matrix::matrix(
	const matrix& A
):
	_lines(A.lines()),
	_columns(A.columns()),
	_block_width(A.block_width()),
	_block_height(A.block_height())
{
	this->_is_transposed = A._is_transposed;
	this->_data = new f32(A.columns()*A.lines());

	memcpy(this->_data, A.data(), sizeof(f32)*A.columns()*A.lines());
}

matrix::matrix(
	const i32 l,
	const i32 c,
	const f32* data,
	const bool transposed,
	const i32 b_w,
	const i32 b_h
): 	_lines(l),
	_columns(c),
	_block_height(b_h), 
	_block_width(b_w), 
	_is_transposed(transposed)
{
	/*
	* MatrixData is stored in a block major fashior for optimizing 
	* cache locality on operations. So the content of the data
	* will be stored as continuously small matrix of size (block_height, block_width)

	*/
	this->_data = new f32[l*c];

	i32 idx = 0;

	for(i32 block_line=0; block_line<l; block_line+=b_h) {
		for(i32 block_column=0; block_column<c; block_column+=b_w) {
			for(i32 y=0; y<b_h; y++) {
				for(i32 x=0; x<b_w; x++) {
					this->_data[idx++] = data[(block_line + y) * c + block_column + x];
				}
			}
		}
	}   
}

f32* matrix::data() const {
	return this->_data;
}

const i32 matrix::lines() const {
	return (1 - this->_is_transposed) * this->_lines + (this->_is_transposed)* this->_columns;
}

const i32 matrix::columns() const {
	return (1 - this->_is_transposed) * this->_columns + (this->_is_transposed)* this->_lines;
}

const i32 matrix::index_of(i32 i, i32 j) const {
	i32 block_line = ((1 - _is_transposed)*i + _is_transposed*j)/_block_height;
	i32 block_column = ((1 - _is_transposed)*j + _is_transposed*i)/_block_width;
	
	i32 y = ((1 - _is_transposed)*i + _is_transposed*j)%_block_height;
	i32 x = ((1 - _is_transposed)*j + _is_transposed*i)%_block_width;

	return block_line*(_columns/_block_width)*(_block_height * _block_width)
		+ block_column *(_block_height * _block_width) + (y*_block_width + x);
}

const i32 matrix::block_width() const {
	return (1 - _is_transposed)*_block_width + _is_transposed*_block_height;
}

const i32 matrix::block_height() const {
	return (1 - _is_transposed)*_block_height + _is_transposed*_block_width;
}

const f32 matrix::get(i32 i, i32 j) const {
	return this->_data[index_of(i,j)];
}

const void matrix::print() {
	for(int i=0; i<lines(); i++) {
		std::cout << "[";
		for(int j=0; j<columns(); j++) {
			std::cout << this->get(i,j);
			if(j != _columns - 1)
				std::cout	<< ", "; 
		}

		std::cout << "]" << std::endl;
	}
	std::cout <<std::endl;
}

matrix matrix::identity(const i32 l, const i32 c, const i32 b_w, const i32 b_h) {
	matrix ident = matrix::zeros(l, c, b_w, b_h);

	for(i32 i=0; i < i32min(l, c); i++)
		ident._data[ident.index_of(i,i)] = 1.f;

	return ident;
}

const void matrix::set(i32 i, i32 j, f32 val) const {
	this->_data[this->index_of(i,j)] = val;
}

void matrix::transpose() {
	this->_is_transposed = true;
}

matrix matrix::operator*(const matrix& B) {
	assert(this->columns() == B.lines());
	assert(this->block_width() == B.block_height());
	
	matrix C = matrix::zeros(this->lines(), B.columns());

	// Iterate over C blocks
	for(i32 i=0; i<this->lines()/this->block_height(); i++) {
		for(i32 j=0; j<B.columns()/B.block_width(); j++) {
			//working C(i,j) block
			// C(i,j) = A(i,:)*B(:,j)
			for(i32 k=0; k<this->columns()/this->block_width(); k++) {
				// k'th A block line
				// k'th B block column

				// Multiply block matrix
				for(i32 y=0; y<B.block_height(); y++) {
					for(i32 x=0; x<this->block_width(); x++) {
						for(i32 q=0; q<this->block_width(); q++) {
							i32 Ai = (B.block_height()*i) % this->lines() + y;
							i32 Aj = (this->block_width()*k) % this->columns() + (x+q) % this->block_width();

							i32 Bi = (B.block_height()*k) % B.lines() + (x+q) % this->block_width();
							i32 Bj = (this->block_width()*j) % B.columns() + x;
							
							i32 Ci = (B.block_height()*i) % this->lines() + y;
							i32 Cj = (this->block_width()*j) % B.columns() + x;
					
							C.set(Ci, Cj, C.get(Ci,Cj) + this->get(Ai, Aj) * B.get(Bi, Bj));
						}
					}
				}
			}
		}
	}	
	return C;
}

matrix matrix::operator+(const matrix& B) {
	assert(this->columns() == B.columns());
	assert(this->lines() == B.lines());

	assert(this->block_height() == B.block_height());
	assert(this->block_width() == B.block_width());

	matrix C = matrix::zeros(this->lines(), B.columns());

	// Iterate over C blocks
	for(i32 i=0; i<this->lines()/this->block_height(); i++) {
		for(i32 j=0; j<B.columns()/B.block_width(); j++) {

			for(i32 x=0; x<B.block_width(); x++) {
				for(i32 y=0; y<this->block_height(); y++) {

					i32 Ci = (this->block_height()*i) % this->lines() + y;
					i32 Cj = (B.block_width()*j) % B.columns() + x;

					C.set(Ci, Cj, this->get(Ci, Cj) + B.get(Ci, Cj));
				}
			}
		}
	}
	return C;
}

matrix matrix::operator-(const matrix& B) {
	assert(this->columns() == B.columns());
	assert(this->lines() == B.lines());

	assert(this->block_height() == B.block_height());
	assert(this->block_width() == B.block_width());

	matrix C = matrix::zeros(this->lines(), B.columns());

	// Iterate over C blocks
	for(i32 i=0; i<this->lines()/this->block_height(); i++) {
		for(i32 j=0; j<B.columns()/B.block_width(); j++) {

			for(i32 x=0; x<B.block_width(); x++) {
				for(i32 y=0; y<this->block_height(); y++) {

					i32 Ci = (this->block_height()*i) % this->lines() + y;
					i32 Cj = (B.block_width()*j) % B.columns() + x;

					C.set(Ci, Cj, this->get(Ci, Cj) - B.get(Ci, Cj));
				}
			}
		}
	}
	return C;
}

matrix matrix::operator/(const matrix& B) {
	assert(this->columns() == B.columns());
	assert(this->lines() == B.lines());

	assert(this->block_height() == B.block_height());
	assert(this->block_width() == B.block_width());

	matrix C = matrix::zeros(this->lines(), B.columns());

	// Iterate over C blocks
	for(i32 i=0; i<this->lines()/this->block_height(); i++) {
		for(i32 j=0; j<B.columns()/B.block_width(); j++) {

			for(i32 x=0; x<B.block_width(); x++) {
				for(i32 y=0; y<this->block_height(); y++) {

					i32 Ci = (this->block_height()*i) % this->lines() + y;
					i32 Cj = (B.block_width()*j) % B.columns() + x;

					C.set(Ci, Cj, this->get(Ci, Cj) / B.get(Ci, Cj));
				}
			}
		}
	}
	return C;
}
