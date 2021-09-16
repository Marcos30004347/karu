#pragma once

#include "algebra/matrix/Matrix.hpp"
#include <initializer_list>

namespace karu::algebra {

class complex 
{
public:
	f32 real;
	f32 imag;

	complex(f32 r, f32 i);

	f32 abs();

	complex conj() const;

	complex operator*(f32 other);
	complex operator+(f32 other);
	complex operator-(f32 other);
	complex operator/(f32 other);

	complex operator*(complex& other);
	complex operator+(complex& other);
	complex operator-(complex& other);
	complex operator/(complex& other);

	complex operator*(const complex& other);
	complex operator+(const complex& other);
	complex operator-(const complex& other);
	complex operator/(const complex& other);
};

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
	Polynomial operator/(const f32 other);
	Polynomial operator+(const Polynomial& other);
	Polynomial operator*(const Polynomial& other);
	Polynomial operator*(const f32 other);
	Polynomial operator-(const Polynomial& other);

	bool operator == (const Polynomial& other);
	
	Polynomial cauchy();
	Polynomial normalized();
	Polynomial dx();
	Polynomial abs();
	Polynomial reverse();
	Polynomial noLeadingTerm();
	Polynomial noConstantTerm();

	f32& operator[](u64 i) const;
	f32 eval(f32 x);
	complex eval(complex x);
	Polynomial root(std::vector<complex>& roots, f32 tol = 0.00001);
	void roots(std::vector<complex>& roots, f32 tol = 0.00001);
};

f32 newtonRaphson(Polynomial& p, f32 tol = 0.0001);

Polynomial addPoly(const Polynomial& A, const Polynomial& B);
Polynomial subPoly(const Polynomial& A, const Polynomial& B);
Polynomial mulPoly(const Polynomial& A, const f32 B);
Polynomial mulPoly(const Polynomial& A, const Polynomial& B);

void divPoly(const Polynomial p,const Polynomial& d, Polynomial& q, Polynomial& r);
void printPoly(Polynomial& p);

bool quadraShift(std::vector<complex>& xs, Polynomial& p, Polynomial& K, complex& x, u64& M, u64& L, bool& qflag, bool& lflag, f32 tol);
bool linearShift(std::vector<complex>& xs, Polynomial& p, Polynomial& K, complex& x, u64 M, u64 L, bool& qflag, bool& lflag, f32 tol);



}
