#pragma once

#include "algebra/matrix/Matrix.hpp"
#include <initializer_list>

namespace karu::algebra {

class Polynomial {
public:
	f32* p_coefficients;
	i64  p_power;

	Polynomial();
	Polynomial(const Polynomial& other);
	Polynomial(u64 power, f32* coefs);
	Polynomial(u64 power, std::initializer_list<f32> coefs);
	~Polynomial();

	i64 power() const;

	Polynomial& operator = (const Polynomial& other);
	Polynomial operator/(const Polynomial& other);
	Polynomial operator+(const Polynomial& other);
	Polynomial operator*(const Polynomial& other);
	Polynomial operator*(const f32 other);
	Polynomial operator-(const Polynomial& other);

	bool operator == (const Polynomial& other);
	
	Polynomial cauchy();
	Polynomial normalized();
	Polynomial dx();
	
	f32& operator[](u64 i) const;
	f32 eval(f32 x);
	f32 root(f32 tol = 0.00001);
	void roots(std::vector<f32>& roots, f32 tol = 0.00001);
};

f32 newtonRaphson(Polynomial& p, f32 tol = 0.0001);

Polynomial addPoly(const Polynomial& A, const Polynomial& B);
Polynomial subPoly(const Polynomial& A, const Polynomial& B);
Polynomial mulPoly(const Polynomial& A, const f32 B);
Polynomial mulPoly(const Polynomial& A, const Polynomial& B);

void divPoly(const Polynomial p,const Polynomial& d, Polynomial& q, Polynomial& r);
void printPoly(Polynomial& p);

}
