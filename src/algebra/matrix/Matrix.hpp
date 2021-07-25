#include "algebra/matrix/MatrixData.hpp"
#include "algebra/vector/Vector.hpp"
#include <initializer_list>

namespace karu::algebra {

class Matrix {
public:
	MatrixData m_data;

	Matrix(const MatrixData& other);
	Matrix(const Vector& other);
	Matrix(const Matrix& other);

	Matrix(u32 l, u32 c);
	Matrix(u32 l, u32 c, f32* data);
	Matrix(u32 l, u32 c, std::initializer_list<f32> data);
	Matrix(u32 l, u32 c, std::initializer_list<f32> data, u32 bx = 16, u32 by = 16);
	Matrix(u32 l, u32 c, u32 bx = 16, u32 by = 16);

	~Matrix();

	const u32 rows() const;
	const u32 columns() const;

	Matrix& operator=(const Vector& other);
	Matrix& operator=(const Matrix& other);
	Matrix operator+(const Matrix& other);
	Matrix operator-(const Matrix& other);
	Matrix operator*(const Matrix& other);

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

}
