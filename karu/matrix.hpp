#ifndef KARU_CORE_MATRIX_H
#define KARU_CORE_MATRIX_H

#include "types.hpp"

namespace karu {

/* 

'block_height' needs to be a divisor of 'lines' and 'block_width' needs to be a
divisor of 'columns', so the matrix can be splited into 
(lines/block_height)*(columns/block_width) submatrices.

Data is stored into a block major fashion, what this means is that 
the first 'block_height*block_width' elements of the array correspont
to a submatrix of size 'block_height*block_width', both blocks and 
submatrices are stored into a row major fashion.
ex:
[1, 2, 3,  4]
[5, 6, 7,  8]
[9,10, 11,12]
[13,14,15,16]
using 2x2 blocks will be stored as:
[1, 2, 5,  6,
	3, 4, 7,  8,
	9,10, 13,14,
	11,12,15,16]

*/

class matrix {
	// storage pointer
	f32* _data;

	// matrix lines
	const i32 _lines;

	// matrix columns
	const i32 _columns;

	// block_height needs to divide columns
	const i32 _block_height;

	// block_width needs to divide lines
	const i32 _block_width;

	// 1 if transposed, 0 otherwise
	bool _is_transposed;

	matrix(
		const i32 lines,
		const i32 columns,
		const i32 block_width = 2,
		const i32 block_height = 2
	);

	const i32 block_width() const;
	const i32 block_height() const;

public:
	matrix(const matrix& A);

	matrix(
		const i32 lines,
		const i32 columns,
		const f32* matrix_data,
		const bool is_transposed = 0,
		const i32 block_width = 2,
		const i32 block_height = 2
	);

	static matrix zeros(
		const i32 lines,
		const i32 columns,
		const i32 block_width = 2,
		const i32 block_height = 2
	);

	static matrix identity(
		const i32 lines,
		const i32 columns,
		const i32 block_width = 2,
		const i32 block_height = 2
	);

	const f32 get(i32 i, i32 j) const;

	const void print();

	const i32 index_of(i32 i, i32 j) const;

	const void set(i32 i, i32 j, f32 val) const;

	f32* data() const;

	const i32 lines() const;

	const i32 columns() const;

	matrix transposed();

	matrix operator*(const matrix& B);
	matrix operator+(const matrix& B);
	matrix operator-(const matrix& B);
	matrix operator/(const matrix& B);
};

};


#endif