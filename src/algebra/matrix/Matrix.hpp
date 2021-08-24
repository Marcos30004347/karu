#pragma once

#include <vector>
#include <initializer_list>
#include "algebra/matrix/MatrixData.hpp"
#include "algebra/vector/Vector.hpp"

namespace karu::algebra {

class Matrix {
public:
	MatrixData m_data;

	Matrix(f32 k);
	Matrix(const MatrixData& other);
	Matrix(const Vector& other);
	Matrix(const Matrix& other);

	Matrix();
	Matrix(u32 l, u32 c, f32* data);
	Matrix(u32 l, u32 c, std::initializer_list<f32> data, u32 bx = 1, u32 by = 1);
	Matrix(u32 l, u32 c, u32 bx = 1, u32 by = 1);

	~Matrix();

	const u32 rows() const;
	const u32 columns() const;

	Matrix& operator=(const Vector& other);
	Matrix& operator=(const Matrix& other);
	Matrix operator+(const Matrix& other);
	Matrix operator-(const Matrix& other);
	Matrix operator*(const Matrix& other);
	Matrix operator*(const f32 other);
	Matrix operator/(const Matrix& other);
	Matrix operator/(const f32 other);

	class MatrixLineGetter {
		friend class Matrix;
		Matrix* parent;
		unsigned int line;
		MatrixLineGetter(Matrix* parent, unsigned int line);
	public:
		f32& operator[](const unsigned int idx);
	};

	MatrixLineGetter operator[](const unsigned int idx);
};

Matrix diag(f32* diag, u32 m, u32 n);
Matrix diag(Matrix& diag, u32 m, u32 n);
Matrix diag(Matrix&& diag, u32 m, u32 n);

Matrix diag(f32* diag, u32 n);
Matrix diag(Matrix& diag);
Matrix diag(Matrix&& diag);

Matrix echelonForm(Matrix matrix);
Matrix nullspace(Matrix matrix);

Matrix transpose(Matrix& matrix);
Matrix transpose(Matrix&& matrix);
Matrix transpose(Matrix* matrix);

Matrix identity(u32 m, u32 n);
void identity(Matrix& I);
void identity(Matrix* I);

std::pair<Matrix, Matrix> LUDecomposition(const Matrix* const A);
std::pair<Matrix, Matrix> LUDecomposition(const Matrix& A);

std::pair<Matrix, Matrix> LUPDecomposition(const Matrix& A);

Matrix LUPSolve(const Matrix& A, const Matrix& P, const Matrix& b);
Matrix LUPInverse(const Matrix& A, const Matrix& P);

f32 LUPDeterminant(const Matrix& A, const Matrix& P);

void printMatrixWithMargin(Matrix& A, unsigned precision = 3, f32 eps = 2.22045e-16);
void printMatrixWithMargin(Matrix* A, unsigned precision = 3, f32 eps = 2.22045e-16);

void printMatrix(Matrix& A, unsigned precision = 5, f32 eps = 2.22045e-16);
void printMatrix(Matrix&& A, unsigned precision = 5, f32 eps = 2.22045e-16);
void printMatrix(Matrix* A, unsigned precision = 5, f32 eps = 2.22045e-16);

}
