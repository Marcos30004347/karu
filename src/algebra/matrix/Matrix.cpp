#include "algebra/matrix/Matrix.hpp"
#include "algebra/matrix/MatrixAdder.hpp"
#include "algebra/matrix/MatrixSubtractor.hpp"
#include "algebra/matrix/MatrixMultiplayer.hpp"

namespace karu::algebra {

Matrix::Matrix(u32 l, u32 c): m_data{l, c, 16, 16} {}
Matrix::Matrix(u32 l, u32 c, f32* data): m_data{l, c, 16, 16, data} {}
Matrix::Matrix(u32 l, u32 c, std::initializer_list<f32> data): m_data{l, c, 16, 16, data} {}
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
	Matrix C(this->m_data.lines(), other.m_data.columns(),  this->m_data.blockWidth(), other.m_data.blockHeight());
	MatrixMultiplayer::mul(&C.m_data, &this->m_data, &other.m_data, false, false);
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

}

