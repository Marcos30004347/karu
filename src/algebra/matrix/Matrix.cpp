#include "algebra/matrix/Matrix.hpp"
#include "algebra/matrix/MatrixAdder.hpp"
#include "algebra/matrix/MatrixSubtractor.hpp"
#include "algebra/matrix/MatrixMultiplayer.hpp"
#include "algebra/matrix/MatrixTransposer.hpp"
#include "algebra/matrix/MatrixDivider.hpp"

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

Matrix transpose(Matrix* other)
{
	Matrix C(other->columns(), other->rows(),  other->m_data.blockHeight(), other->m_data.blockWidth());
	MatrixTransposer::transpose(&C.m_data, &other->m_data);
	return C;
}

void printMatrixWithMargin(Matrix* A)
{
	for(int i=0;i<A->m_data.m_stored_lines; i++)
	{
		for(int j=0; j<A->m_data.m_stored_column;j++)
		{
			std::cout << A->m_data.get(i,j) << " ";
		}
		std::cout << std::endl;
	}
}

void printMatrixWithMargin(Matrix& A)
{
	for(int i=0;i<A.m_data.m_stored_lines; i++)
	{
		for(int j=0; j<A.m_data.m_stored_column;j++)
		{
			std::cout << A.m_data.get(i,j) << " ";
		}
		std::cout << std::endl;
	}
}

void printMatrix(Matrix* A)
{
	for(int i=0;i<A->m_data.m_lines; i++)
	{
		for(int j=0; j<A->m_data.m_columns;j++)
		{
			std::cout << A->m_data.get(i,j) << " ";
		}
		std::cout << std::endl;
	}
}

void printMatrix(Matrix& A)
{
	for(int i=0;i<A.m_data.m_lines; i++)
	{
		for(int j=0; j<A.m_data.m_columns;j++)
		{
			std::cout << A.m_data.get(i,j) << " ";
		}
		std::cout << std::endl;
	}
}
}

