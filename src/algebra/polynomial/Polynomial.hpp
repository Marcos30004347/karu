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
	complex(f32 r);
	complex();

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

class RealPoly {
public:
	f32* p_coefficients;
	i64  p_power;

	RealPoly();
	RealPoly(const RealPoly& other);
	RealPoly(u64 power);
	RealPoly(u64 power, f32* coefs);
	RealPoly(u64 power, std::initializer_list<f32> coefs);
	~RealPoly();

	i64 power() const;

	RealPoly& operator = (const RealPoly& other);
	RealPoly operator/(const RealPoly& other);
	RealPoly operator/(const f32 other);
	RealPoly operator+(const RealPoly& other);
	RealPoly operator*(const RealPoly& other);
	RealPoly operator*(const f32 other);
	RealPoly operator-(const RealPoly& other);

	bool operator == (const RealPoly& other);
	
	RealPoly cauchy();
	RealPoly normalized();
	RealPoly dx();
	RealPoly abs();
	RealPoly reverse();
	RealPoly noLeadingTerm();
	RealPoly noConstantTerm();

	f32& operator[](u64 i) const;
	f32 eval(f32 x);
	complex eval(complex x);
	RealPoly root(std::vector<complex>& roots, f32 tol = 0.00001);
	// void roots(std::vector<complex>& roots, f32 tol = 0.00001);
};

f32 newtonRaphson(RealPoly& p, f32 tol = 0.0001);

RealPoly addPoly(const RealPoly& A, const RealPoly& B);
RealPoly subPoly(const RealPoly& A, const RealPoly& B);
RealPoly mulPoly(const RealPoly& A, const f32 B);
RealPoly mulPoly(const RealPoly& A, const RealPoly& B);
std::vector<complex> polyRoots(RealPoly p, f32 tol = 2.22e-16);

void divPoly(const RealPoly p,const RealPoly& d, RealPoly& q, RealPoly& r);
void printPoly(RealPoly& p);

bool quadraShift(std::vector<complex>& xs, RealPoly& p, RealPoly& K, complex& x, u64& M, u64& L, bool& qflag, bool& lflag, f32 tol);
bool linearShift(std::vector<complex>& xs, RealPoly& p, RealPoly& K, complex& x, u64 M, u64 L, bool& qflag, bool& lflag, f32 tol);



}
