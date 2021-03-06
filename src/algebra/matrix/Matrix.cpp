#include "algebra/matrix/Matrix.hpp"
#include "algebra/matrix/MatrixAdder.hpp"
#include "algebra/matrix/MatrixSubtractor.hpp"
#include "algebra/matrix/MatrixMultiplayer.hpp"
#include "algebra/matrix/MatrixTransposer.hpp"
#include "algebra/matrix/MatrixDivider.hpp"
#include "algebra/matrix/MatrixLU.hpp"
#include "algebra/matrix/MatrixEchelonForm.hpp"
#include "algebra/matrix/MatrixNullSpace.hpp"
#include "algebra/linear/Linear.hpp"
#include <assert.h>
#include <cmath>
#include <iomanip>

namespace karu::algebra {


Matrix::Matrix(f32 k): m_data{1, 1, 1, 1, {k}} {}
Matrix::Matrix(): m_data{0, 0, 0, 0} {}
Matrix::Matrix(u32 l, u32 c, f32* data): m_data{l, c, 1, 1, data} {}
Matrix::Matrix(u32 l, u32 c, std::initializer_list<f32> data, u32 bx, u32 by): m_data{l, c, bx, by, data} {}
Matrix::Matrix(u32 l, u32 c, u32 bx, u32 by): m_data{l, c, bx, by} {}

Matrix::Matrix(const Vector& other)
{
	u32 size = other.size();
	this->m_data = MatrixData(size, 1, 1, 1);
}

Matrix::Matrix(const MatrixData& other)
{
	this->m_data = other;
}

Matrix::Matrix(const Matrix& other)
{
	this->m_data = other.m_data;
}

Matrix& Matrix::operator=(const Vector& other)
{

	u32 size = other.size();
	this->m_data = MatrixData(size, 1, 1, 1);
	return *this;
}

Matrix& Matrix::operator=(const Matrix& other)
{

	this->m_data = other.m_data;
	return *this;
}

Matrix Matrix::operator+(const Matrix& other)
{
	Matrix C(this->m_data.lines(), this->m_data.columns(),  this->m_data.blockHeight(), other.m_data.blockWidth());
	MatrixAdder::add(&C.m_data, &this->m_data, &other.m_data, false, false);
	return C;
}

Matrix Matrix::operator-(const Matrix& other)
{
	Matrix C(this->m_data.lines(), this->m_data.columns(),  this->m_data.blockHeight(), other.m_data.blockWidth());
	MatrixSubtractor::sub(&C.m_data, &this->m_data, &other.m_data, false, false);
	return C;
}

Matrix Matrix::operator*(const Matrix& other)
{
	// std::cout << "ASDASD\n";
	if(other.columns() == 1 && other.rows() == 1)
	{
		Matrix D(this->rows(), this->columns(),  this->m_data.blockWidth(), this->m_data.blockHeight());
		MatrixMultiplayer::mul(&D.m_data, &this->m_data, other.m_data.get(0,0));
		return D;
	}

	if(this->columns() == 1 && this->rows() == 1)
	{
		Matrix D(other.rows(), other.columns(),  other.m_data.blockWidth(), other.m_data.blockHeight());
		MatrixMultiplayer::mul(&D.m_data, &other.m_data, this->m_data.get(0,0));
		return D;
	}

	Matrix C(this->m_data.lines(), other.m_data.columns(),  this->m_data.blockWidth(), other.m_data.blockHeight());
	MatrixMultiplayer::mul(&C.m_data, &this->m_data, &other.m_data, false, false);
	return C;
}

Matrix Matrix::operator/(const Matrix& other)
{
	if(other.columns() == 1 && other.rows() == 1)
	{
		Matrix D(this->rows(), this->columns(),  this->m_data.blockWidth(), this->m_data.blockHeight());
		MatrixDivider::div(&D.m_data, &this->m_data, other.m_data.get(0,0));
		return D;
	}

	if(this->columns() == 1 && this->rows() == 1)
	{
		Matrix D(other.rows(), other.columns(),  other.m_data.blockWidth(), other.m_data.blockHeight());
		MatrixDivider::div(&D.m_data, &other.m_data, this->m_data.get(0,0));
		return D;
	}

	Matrix C(this->m_data.lines(), other.m_data.columns(),  this->m_data.blockWidth(), other.m_data.blockHeight());
	MatrixDivider::div(&C.m_data, &this->m_data, &other.m_data, false, false);
	return C;
}

Matrix Matrix::operator/(const f32 other)
{
	Matrix C(this->m_data.lines(), this->m_data.columns(),  this->m_data.blockWidth(), this->m_data.blockHeight());
	MatrixDivider::div(&C.m_data, &this->m_data, other);
	return C;
}

Matrix Matrix::operator*(const f32 other)
{
	Matrix C(this->m_data.lines(), this->m_data.columns(),  this->m_data.blockWidth(), this->m_data.blockHeight());
	MatrixMultiplayer::mul(&C.m_data, &this->m_data, other);
	return C;
}


Matrix::MatrixLineGetter::MatrixLineGetter(Matrix* parent, unsigned int line)
{
	this->line = line;
	this->parent = parent;
}

f32& Matrix::MatrixLineGetter::operator[](const unsigned int idx)
{
	return parent->m_data.get(this->line, idx);
}

Matrix::MatrixLineGetter Matrix::operator[](const unsigned int idx)
{
	return Matrix::MatrixLineGetter(this, idx);
}

Matrix::~Matrix(){}

const u32 Matrix::rows() const
{
	return this->m_data.lines();
}

const u32 Matrix::columns() const
{
	return this->m_data.columns();
}

Matrix transpose(Matrix& other)
{
	Matrix C(other.columns(), other.rows(),  other.m_data.blockHeight(), other.m_data.blockWidth());
	MatrixTransposer::transpose(&C.m_data, &other.m_data);
	return C;
}

Matrix transpose(Matrix&& other)
{
	Matrix C(other.columns(), other.rows(),  other.m_data.blockHeight(), other.m_data.blockWidth());
	MatrixTransposer::transpose(&C.m_data, &other.m_data);
	return C;
}

Matrix transpose(Matrix* other)
{
	Matrix C(other->columns(), other->rows(),  other->m_data.blockHeight(), other->m_data.blockWidth());
	MatrixTransposer::transpose(&C.m_data, &other->m_data);
	return C;
}

std::pair<Matrix, Matrix> LUDecomposition(const Matrix* const A)
{
	Matrix L(A->rows(), A->columns(), A->m_data.blockWidth(), A->m_data.blockHeight());
	Matrix U(A->rows(), A->columns(), A->m_data.blockWidth(), A->m_data.blockHeight());
	MatrixLU::LUdecompose(&L.m_data, &U.m_data, &A->m_data);
	return {L, U};
}

std::pair<Matrix, Matrix> LUDecomposition(const Matrix& A)
{
	Matrix L(A.rows(), A.columns(), A.m_data.blockWidth(), A.m_data.blockHeight());
	Matrix U(A.rows(), A.columns(), A.m_data.blockWidth(), A.m_data.blockHeight());
	MatrixLU::LUdecompose(&L.m_data, &U.m_data, &A.m_data);
	return {L, U};
}

std::pair<Matrix, Matrix> LUPDecomposition(const Matrix& A)
{
	Matrix R = Matrix(A);
	Matrix P(A.rows()+1, 1);
	MatrixLU::LUPdecompose(&R.m_data, &P.m_data);
	return {R, P};
}

Matrix LUPSolve(const Matrix& A, const Matrix& P, const Matrix& b)
{
	Matrix x(A.rows(), 1);
	MatrixLU::LUPSolve(&A.m_data, &P.m_data, &b.m_data, &x.m_data);
	return x;
}

Matrix LUPInverse(const Matrix& A, const Matrix& P)
{
	Matrix Inv(A.rows(), A.columns(), A.m_data.blockWidth(), A.m_data.blockHeight());
	MatrixLU::LUPInvet(&A.m_data, &P.m_data, &Inv.m_data);
	return Inv;
}

f32 LUPDeterminant(const Matrix& A, const Matrix& P)
{
	return MatrixLU::LUPDeterminant(&A.m_data, &P.m_data);
}

Matrix echelonForm(Matrix matrix)
{
	MatrixEchelonForm::toEchelonForm(&matrix.m_data);
	return matrix;
}

Matrix nullspace(Matrix matrix)
{
	MatrixData null;
	MatrixNullSpace::nullSpace(&matrix.m_data, null);
	return Matrix(null);
}

Matrix diag(Matrix&& diag) 
{
	Matrix a(diag.rows(),diag.rows());
	for(i32 i=0; i<diag.rows(); i++)
		a[i][i] = diag[i][0];
	return a;
}

Matrix diag(Matrix& diag) 
{
	Matrix a(diag.rows(),diag.rows());
	for(i32 i=0; i<diag.rows(); i++)
		a[i][i] = diag[i][0];
	return a;
}

Matrix diag(f32* diag, u32 n) 
{
	Matrix a(n,n);
	for(i32 i=0; i<n; i++)
		a[i][i] = diag[i];
	return a;
}

Matrix diag(f32* diag, u32 m, u32 n) 
{
	Matrix a(m,n);

	for(i32 i=0; i<std::min(n, m); i++)
		a[i][i] = diag[i];

	return a;
}

Matrix bidiag(f32* diag, f32* sdiag, u32 m, u32 n) 
{
	Matrix a(m,n);

	for(i32 i=0; i<std::min(n, m); i++)
	{
		a[i][i] = diag[i];
	
		if(i < std::min(n,m) - 1)
		{
			a[i][i+1] = sdiag[i+1];
		}
	}

	return a;
}

Matrix diag(Matrix&& diag, u32 m, u32 n) 
{
	assert(diag.rows() == n);
	Matrix a(m,n);
	for(i32 i=0; i<std::min(n, m); i++)
		a[i][i] = diag[i][0];
	return a;
}

Matrix diag(Matrix& diag, u32 m, u32 n) 
{
	assert(diag.rows() == n);

	Matrix a(m,n);
	for(i32 i=0; i<n; i++)
		a[i][i] = diag[i][0];
	return a;
}

void printMatrixWithMargin(Matrix* A, unsigned precision, f32 eps)
{
	for(int i=0;i<A->m_data.m_stored_lines; i++)
	{
		for(int j=0; j<A->m_data.m_stored_column;j++)
		{
			f32 val = fabs(A->m_data.get(i,j)) > eps ? A->m_data.get(i,j) : 0.0;
			std::cout << std::fixed << std::setprecision(precision) << std::setw(precision+8) << val << " ";
		}
		std::cout << std::endl;
	}
}

void printMatrixWithMargin(Matrix& A, unsigned precision, f32 eps)
{
	for(int i=0;i<A.m_data.m_stored_lines; i++)
	{
		for(int j=0; j<A.m_data.m_stored_column;j++)
		{
			f32 val = fabs(A.m_data.get(i,j)) > eps ? A.m_data.get(i,j) : 0.0;

			std::cout << std::fixed << std::setprecision(precision) << std::setw(precision+8) <<val << " ";
		}
		std::cout << std::endl;
	}
}

void printMatrix(Matrix* A, unsigned precision, f32 eps)
{
	for(int i=0;i<A->m_data.m_lines; i++)
	{
		for(int j=0; j<A->m_data.m_columns;j++)
		{
			f32 val = fabs(A->m_data.get(i,j)) > eps ? A->m_data.get(i,j) : 0.0;

			std::cout << std::scientific << std::setprecision(precision) << std::setw(precision+8) << val << ", ";
		}
		std::cout << std::endl;
	}
}

void printMatrix(Matrix& A, unsigned precision, f32 eps)
{
	for(int i=0;i<A.m_data.m_lines; i++)
	{
		for(int j=0; j<A.m_data.m_columns;j++)
		{
			f32 val = fabs(A.m_data.get(i,j)) > eps ? A.m_data.get(i,j) : 0.0;

			std::cout << std::scientific << std::setprecision(precision) << std::setw(precision+8) <<val << ", ";
		}
		std::cout << std::endl;
	}
}

void printSubMatrix(Matrix& A, i32 p, i32 q, i32 r, i32 s, unsigned precision, f32 eps)
{
	for(int i=p;i<q; i++)
	{
		for(int j=r; j<s;j++)
		{
			f32 val = fabs(A.m_data.get(i,j)) > eps ? A.m_data.get(i,j) : 0.0;

			std::cout << std::scientific << std::setprecision(precision) << std::setw(precision+8) <<val << ", ";
		}
		std::cout << std::endl;
	}
}

void printMatrix(Matrix&& A, unsigned precision, f32 eps)
{
	for(int i=0;i<A.m_data.m_lines; i++)
	{
		for(int j=0; j<A.m_data.m_columns;j++)
		{			
			f32 val = fabs(A.m_data.get(i,j)) > eps ? A.m_data.get(i,j) : 0.0;
		
			std::cout << std::scientific << std::setprecision(precision) << std::setw(precision+8) <<val << ", ";
		}
		std::cout << std::endl;
	}
}

Matrix identity(u32 m, u32 n)
{
	Matrix I(m,n);

	for(i32 i = 0; i < std::min(n, m); i++)
	{
		I[i][i] = 1.0;
	}

	return I;
}

}

