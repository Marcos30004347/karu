#include "algebra/polynomial/Polynomial.hpp"
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>

namespace karu::algebra {


complex::complex(f32 r, f32 i)
{
	this->real = r;
	this->imag = i;
}

complex complex::conj() const
{
	return complex(this->real, -1 * this->imag);
}

complex complex::operator*(complex& other)
{
	// 1 = i^2
	return complex(
		(this->real * other.real) + (-1 * this->imag * other.imag),
		(this->real * other.imag) + (this->imag * other.real)
	);
}

complex complex::operator+(complex& other)
{
	return complex(
		this->real + other.real,
		this->imag + other.imag
	);
}

complex complex::operator-(complex& other)
{
	return complex(
		this->real - other.real,
		this->imag - other.imag
	);
}

complex complex::operator/(complex& other)
{
	complex a = *this;

	complex b = other;

	complex bc = other.conj();

	complex n = a * bc;
	complex d = b * bc;

	return complex(n.real / d.real, n.imag / d.real);
}



complex complex::operator*(const complex& other)
{
	// 1 = i^2
	return complex(
		(this->real * other.real) + (-1 * this->imag * other.imag),
		(this->real * other.imag) + (this->imag * other.real)
	);
}
complex complex::operator+(const complex& other)
{
	return complex(
		this->real + other.real,
		this->imag + other.imag
	);
}

complex complex::operator-(const complex& other)
{
	return complex(
		this->real - other.real,
		this->imag - other.imag
	);
}

complex complex::operator/(const complex& other)
{
	complex a = *this;

	complex b = other;

	complex bc = other.conj();

	complex n = a * other.conj();
	complex d = b * other.conj();

	return complex(n.real / d.real, n.imag / d.real);
}


f32 complex::abs()
{
	return std::hypot(this->real, this->imag);
}
complex complex::operator*(f32 other)
{
	return complex(this->real * other, this->imag * other);
}

complex complex::operator+(f32 other)
{
	return complex(this->real + other, this->imag + other);
}

complex complex::operator-(f32 other)
{
	return complex(this->real - other, this->imag - other);
}

complex complex::operator/(f32 other)
{
	return complex(this->real / other, this->imag / other);
}

void printPoly(Polynomial& p)
{

	for(i64 i=0;i<p.power()+1; i++)
	{
		if(i > 1)
			printf("%g*x^%lli ", p[i], i);
		else if(i == 1)
			printf("%g*x ", p[i]);
		else
			printf("%g ", p[i]);
		if(i!=p.power())
			printf("+ ");
	}

	std::cout << "\n";
}

void divPoly(const Polynomial p, const Polynomial& d, Polynomial& q, Polynomial& r)
{
	i64 i,j;
	i64 power = p.power() - d.power();
	f32 ratio;

	if(power < 0) return;

	q = Polynomial(power, {});

	r = Polynomial(p);

	for(i=p.power(); i>=d.power(); i--)
	{
		q[i - d.power()] = ratio = r[i]/d[d.p_power];

		// std::cout << "quotient: " << r[i]<< std::endl;
		// std::cout << "quotient: " << d[d.p_power] << std::endl;

		r[i] = 0;

		for(j=0; j<d.power(); j++)
		{
			r[i - d.power() + j] -= d[j]*ratio;
		}
	}
	while(r.p_power >= 0 && r[--r.p_power] == 0);
}


void divPoly(const Polynomial p, const f32 d, Polynomial& q)
{
	q = Polynomial(p);
	for (int j=0; j<=p.power(); j++)
		q[j] = p[j]/d;
}

Polynomial addPoly(const Polynomial& A, const Polynomial& B)
{
	Polynomial X = A.power() > B.power() ? Polynomial(A) : Polynomial(B);
	Polynomial Y = A.power() > B.power() ? Polynomial(B) : Polynomial(A);

	for(i64 i=0; i<=Y.power(); i++)
	{
		X[i] += Y[i];
	}

	while(X[X.p_power] == 0) { X.p_power--; }
	return X;
}

Polynomial subPoly(const Polynomial& A, const Polynomial& B)
{
	Polynomial B_ = mulPoly(B, -1.0);

	Polynomial X = A.power() > B.power() ? Polynomial(A) : Polynomial(B_);
	Polynomial Y = A.power() > B.power() ? Polynomial(B_) : Polynomial(A);

	for(i64 i=0; i<=Y.power(); i++)
	{
		X[i] += Y[i];
	}

	while(X[X.p_power] == 0) { X.p_power--; }

	return X;
}

Polynomial mulPoly(const Polynomial& A, const Polynomial& B)
{
	u64 m = A.power() + 1;
	u64 n = B.power() + 1;

	Polynomial AB(m+n-2,{});

	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
		{
			AB[i+j] += A[i]*B[j];
		}
	}
	return AB;
}

Polynomial mulPoly(const Polynomial& A, f32 B)
{
	Polynomial AB(A);
	for (int j=0; j<=AB.power(); j++)
		AB[j] = A[j]*B;
	return AB;
}

Polynomial::Polynomial()
{
	this->p_power = 0;
	this->p_coefficients = new f32[this->p_power+1];
	this->p_coefficients[0] = 0.0;
}

Polynomial::Polynomial(const Polynomial& other)
{
	this->p_power = other.p_power;
	this->p_coefficients = new f32[this->p_power+1];
	std::copy(other.p_coefficients, other.p_coefficients + other.p_power+1, this->p_coefficients);
}

Polynomial::Polynomial(u64 power, std::initializer_list<f32> coefs)
{
	i64 i;
	this->p_power = power;
	this->p_coefficients = new f32[++power];

	if(coefs.size() == 0)
	{
		std::fill(&this->p_coefficients[0], &this->p_coefficients[power-1]+1, 0);
	}
	else
	{
		std::copy(coefs.begin(), coefs.end(), this->p_coefficients);
	}
}

Polynomial::Polynomial(u64 power, f32* coefs)
{
	i64 i;

	this->p_power = power;
	this->p_coefficients = new f32[++power];

	if(coefs == nullptr)
	{
		std::fill(&this->p_coefficients[0], &this->p_coefficients[power-1]+1, 0);
	}
	else
	{
		std::copy(&coefs[0], &coefs[power-1]+1, this->p_coefficients);
	}
}

Polynomial::~Polynomial(){
	delete this->p_coefficients;
}

i64 Polynomial::power() const
{
	return this->p_power;
}


Polynomial Polynomial::operator/(const Polynomial& other)
{
	Polynomial q,r;
	divPoly(*this, other, q, r);
	return q;
}

Polynomial Polynomial::operator/(const f32 other)
{
	Polynomial q;
	divPoly(*this, other, q);
	return q;
}

Polynomial Polynomial::operator*(const f32 other)
{
	return mulPoly(*this, other);
}

Polynomial Polynomial::operator*(const Polynomial& other)
{
	return mulPoly(*this, other);
}

Polynomial Polynomial::operator+(const Polynomial& other)
{
	return addPoly(*this, other);
}

Polynomial Polynomial::operator-(const Polynomial& other)
{
	return subPoly(*this, other);
}

Polynomial& Polynomial::operator = (const Polynomial& other)
{
	this->p_power = other.p_power;
	this->p_coefficients = new f32[this->p_power+1];
	std::copy(other.p_coefficients, other.p_coefficients + other.p_power+1, this->p_coefficients);
	return *this;
}

bool Polynomial::operator == (const Polynomial& other)
{
	if(other.p_power != this->p_power)
		return false;

	for(i64 i=0; i<=this->p_power; i++)
	{
		if(this->p_coefficients[i] != other.p_coefficients[i])
			return false;
	}

	return true;
}

f32& Polynomial::operator[](u64 i) const
{
	return this->p_coefficients[i];
}

// horners method
f32 Polynomial::eval(f32 x)
{
	f32 px = this->p_coefficients[this->power()];

	for(i64 i = this->power() - 1; i >= 0; i--)
	{
		px = px * x + this->p_coefficients[i];
	}

	return px;
}

complex pow(complex c, u32 i)
{
	complex p = c;

	for(i32 j=0; j<i; j++)
	{
		p = p * c;
	}

	return p;
}

// horners method
complex Polynomial::eval(complex x)
{
	complex px(0, 0);

	for(i64 i = this->power(); i >= 0; i--)
	{
		px = px * x + this->p_coefficients[i];
	}

	return px;
}

Polynomial Polynomial::dx()
{
	Polynomial d(this->power() - 1, {});

	for(i64 i=1; i<=this->power(); i++)
	{
		d[i-1] = this->p_coefficients[i]*i;
	}

	return d;
}

Polynomial Polynomial::abs()
{
	Polynomial d(this->power(), {});

	for(i64 i=0; i<=this->power(); i++)
	{
		d[i] = this->p_coefficients[i] >= 0 ? this->p_coefficients[i] : -this->p_coefficients[i];
	}

	return d;
}

Polynomial Polynomial::reverse()
{
	Polynomial d(this->power(), {});

	for(i64 i=0; i<=this->power(); i++)
	{
		d[this->power() - i] = this->p_coefficients[i];
	}

	return d;
}

Polynomial Polynomial::noLeadingTerm()
{
	Polynomial d(this->power() - 1, {});

	for(i64 i=0; i<this->power(); i++)
	{
		d[i] = this->p_coefficients[i];
	}


	return d;
}

Polynomial Polynomial::noConstantTerm()
{
	Polynomial d(this->power(), {});

	for(i64 i=0; i<=this->power(); i++)
	{
		d[i] = this->p_coefficients[i];
	}

	d[0] = 0;

	return d;
}

Polynomial Polynomial::normalized()
{

	Polynomial n(this->power(), {});

	if(n[n.power()] == 1)
		return n;

	for(i64 i=0; i<=this->power(); i++)
	{
		n[i] = this->p_coefficients[i]/this->p_coefficients[this->p_power];
	}

	return n;
}

Polynomial Polynomial::cauchy()
{
	Polynomial p = this->normalized();

	for(i64 i=0; i<=p.power(); i++)
	{
		p[i] = fabs(p[i]);
	}

	p[p.power()] *= -1;

	return p;
}

f32 randomBetweenZeroOne()
{
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	return unif(rng);
}

f32 newtonRaphson(Polynomial& p, f32 x0, f32 tol, u32 max_it = 1000000)
{
	Polynomial pdx = p.dx();

	f32 x = x0;

	while(pdx.eval(x) == 0)
	{
		x += 0.01;//randomBetweenZeroOne();
	}
	f32 tmp;

	for(i32 i =0; i<max_it; i++)
	{
		tmp = x;

		x = x - p.eval(x)/pdx.eval(x);

		if(fabs(p.eval(x)) <= tol)
		{
			break;
		}
	}

	return x;
}

f32 rootRadius(Polynomial& p, f32 tol = 0.001)
{
	Polynomial p_ = p.abs();

	p_[0] *= -1;

	return newtonRaphson(p_, 1.0, tol, 100);
}

void jenkinsTraubPhaseOne(std::vector<complex>& roots, Polynomial& p, Polynomial& K, complex& s, u64 M, f32 tol)
{
	Polynomial z = Polynomial(1, { 0,1 });

	K = p.dx() / (p.power() + 1);

	for(i64 i=1; i<M; i++)
	{
		K = (K.noConstantTerm() + p.noConstantTerm()*(-K[0]/p[0]))/z;
	}
}

Polynomial nextSigma(Polynomial& p, Polynomial& sigma, Polynomial& K, f32& a, f32& b, f32& c, f32& d)
{
	f32 u = sigma[sigma.power() - 1];
	f32 v = sigma[sigma.power() - 2];

	f32 b1 = -K[0] / p[0];
	f32 b2 = -K[1] + b1 * p[1] / p[0];

	f32 a1 = b * c - a * d;
	f32 a2 = a * c + u * a * d + v * b * d;
	f32 c2 = b1 * a2;
	f32 c3 = b1 * b1 * (a * a + u * a * b + v * b * b);
	f32 c4 = v * b2 * a1 - c2 - c3;
	f32 c1 = c * c + u * c * d + v * d * d + b1 * (a * c + u * b * c + v * b * d) - c4;
	f32 du = -(u * (c2 + c3) + v * (b1 * a1 + b2 * a2)) / c1;
	f32 dv = v * c4 / c1;

	return Polynomial(2, {v + dv, u + du, 1});
}

Polynomial nextKQuaraticShift(Polynomial K, Polynomial& sigma, Polynomial& p_q, Polynomial& k_q, f32& a, f32& b, f32& c, f32& d)
{
	f32 t = (a * a + sigma[sigma.power() - 1] * a * b + sigma[sigma.power() - 2] * b * b) / (b * c - a * d);

	Polynomial lin(2, {-(a * c + sigma[sigma.power() - 1] * a * d + sigma[sigma.power() - 2] * b * d) / (b * c - a * d), 1});

	K = k_q * t + lin * p_q;

	K[0] += b;

	return K;
}

int fixedShiftKPoly(Polynomial& p, Polynomial& sigma, Polynomial& K, complex& x, f32& a, f32& b, f32& c, f32& d, u32 max_it)
{
	Polynomial p_r, p_q, k_q, k_r;

	sigma[0] = x.real*x.real + x.imag*x.imag;
	sigma[1] = -2.0 * x.real;
	sigma[2] = 1;

	divPoly(p, sigma, p_q, p_r);

	b = p_r[p_r.power()];
	a = p_r[p_r.power() - 1] - b * sigma[sigma.power() - 1];

	complex p_x = x.conj() * -b + a;

	complex t_l[3] = {complex(0,0),complex(0,0),complex(0,0)};
	f32 s_l[3] = {0,0,0};

	for(i32 j=0; j < max_it; j++)
	{
		K = K.normalized();

		divPoly(K, sigma, k_q, k_r);

		d = k_r[k_r.power()];
		c = k_r[k_r.power() - 1] - d - sigma[sigma.power() - 1];

		Polynomial s = nextSigma(p, sigma, K, a, b, c, d);

		complex k_x = x.conj() * -d + c;

		t_l[0] = t_l[1];
		s_l[0] = s_l[1];

		t_l[1] = t_l[2];
		s_l[1] = s_l[2];

		t_l[2] = x - p_x / k_x;
		s_l[2] = s[2];

		// return if converge
		if(
			std::fabs(s_l[1] - s_l[0]) < std::fabs(s_l[0]) / 2.0 &&
			std::fabs(s_l[2] - s_l[1]) < std::fabs(s_l[1]) / 2.0
		) return 1;
		if(
			(t_l[1] - t_l[0]).abs() < (t_l[0]).abs() / 2.0 &&
			(t_l[2] - t_l[1]).abs() < (t_l[1]).abs() / 2.0
		) return 2;

		K = nextKQuaraticShift(K, sigma, p_q, k_q, a, b, c, d);
	}

	return 0;
}


int jenkinsTraubPhaseTwo(std::vector<complex>& roots, Polynomial& p, Polynomial& K, complex& x, u64 L, f32 tol)
{
	u32 i, j;

	f32 a, b, d, c;

	f32 deg2rad = M_PI / 180.0;
  f32 phi = 49.0 * deg2rad;
	f32 radius = rootRadius(p, 1e-2);
	std::cout << "root_radius" << "\n";
	std::cout << radius << "\n";
	Polynomial sigma(2, {});

	for(i=0; i<200; i++)
	{
		x = complex(radius, 0) * complex(cos(phi), sin(phi));
	
		std::cout << "root" << "\n";
		std::cout << x.real << " + " << x.imag << "i\n";
	
		i32 j = fixedShiftKPoly(p, sigma, K, x, a, b, c, d, L * (i + 1));
		
		std::cout << j << "\n";
		std::cout << "phase2 K:\n";
		printPoly(K);
		
		if(j != 0)
		{
			return j;
		}

		phi += 94.0 * deg2rad;
	}

	return 0;
}

bool hasConverged(std::vector<complex>& roots, f32 tol)
{
	if(roots.size() != 3) return false;

	f32 e[2];

	e[1] = (roots[2] - roots[1]).abs();
	e[0] = (roots[1] - roots[0]).abs();

	f32 r = roots[1].abs();

	if(e[1] <= e[0])
	{
		if(r <= tol)
		{
			return e[1] < tol;
		}
		else
		{
			return e[1] / r <= tol;
		}
	}

	return false;
}

bool hasConverged(std::vector<f32>& roots, f32 tol)
{
	if(roots.size() != 3) return false;

	f32 e[2];

	e[1] = fabs(roots[2] - roots[1]);
	e[0] = fabs(roots[1] - roots[0]);

	f32 r = fabs(roots[1]);

	if(e[1] <= e[0])
	{
		if(r <= tol)
		{
			return e[1] < tol;
		}
		else
		{
			return e[1] / r <= tol;
		}
	}

	return false;
}
void findQuadraticRoots(Polynomial& p, complex& x1, complex& x2)
{
	f32 a = p[2], b = p[1], c = p[0];

	f32 d = b * b - 4 * a * c;

	if(d >= 0)
	{
		if(b >= 0)
		{
			x1 = complex((-b - sqrt(fabs(d))) / (2.0 * a), 0);
			x2 = complex((2.0 * c) / (-b - sqrt(fabs(d))), 0);
		}
		else
		{
			x1 = complex((2.0 * c) / (-b + sqrt(fabs(d))), 0);
			x2 = complex((-b + sqrt(fabs(d))) / (2.0 * a), 0);
		}
	}
	else
	{
		x1 = complex(-b / (2.0 * a), sqrt(fabs(d)) / (2.0 * a));
		x2 = complex(-b / (2.0 * a), -sqrt(fabs(d)) / (2.0 * a));
	}

}


bool linearShift(std::vector<complex>& zeros, Polynomial& p, Polynomial& K, complex& x, u64 M, u64 L, bool& qflag, bool& lflag, f32 tol)
{
	if(lflag)
	{
		return false;
	}

	std::vector<f32> roots;

	Polynomial defl_p, defl_k;
	
	f32 p_x = 0, px = 0, kx = 0, d = 0;

	// s^(L) = Re(s1 - P(s1)/K^(L)(s1))
	f32 s = (x - p.eval(x) / K.eval(x)).real;

	roots.push_back(s);

	for(i32 i = 0; i < L; i++)
	{
		if(hasConverged(roots, tol))
		{
			// std::cout << "1#\n";
			zeros.push_back(complex(roots[1], 0));
			p = defl_p;
			return true;
		}

		p_x = px;

		defl_p = p / Polynomial(1, {-s, 1});
		px = p.eval(s);

		if(fabs(px) <= tol)
		{
			// std::cout << "2#\n";
	
			zeros.push_back(complex(s, 0));
			p = defl_p;
	
			return true;
		}
	
		defl_k = K / Polynomial(1, { -s, 1 });

		kx = K.eval(s);
		
		K = defl_k +  defl_p * -kx / px;

		K = K.normalized();

		kx = K.eval(s);
	
		d = px / kx;
	
		s = s - px / kx;

		roots.push_back(s);

		if(roots.size() > 3)
		{
			roots.erase(roots.begin());
		}

		if(i >= 2 && fabs(d) < 0.001 * fabs(s) && fabs(p_x) < fabs(p_x))
		{
			// std::cout << "4#\n";
			complex root(s, 0);
			return quadraShift(zeros, p, K, root, M, L, qflag, lflag, tol);
		}
	}

	lflag = true;
	// std::cout << "5#\n";
	return quadraShift(zeros, p, K, x, M, L, qflag, lflag, tol);
}

bool quadraShift(std::vector<complex>& zeros, Polynomial& p, Polynomial& K, complex& x, u64& M, u64& L, bool& qflag, bool& lflag, f32 tol)
{
	if(qflag)
	{
		return false;
	}

	f32 step = 0.01;
	f32 a, b, c, d;
	f32 px = 0, p_x = 0, tmp = 0;

	std::vector<complex> roots[2];
	
	Polynomial p_q, p_r, k_q, k_r, sigma(2, {});

	sigma[0] = x.real * x.real + x.imag * x.imag;
	sigma[1] = -2.0 * x.real;
	sigma[2] = 1;

	roots[0].push_back(x);
	roots[1].push_back(x.conj());

	bool fixed_shift = false;

	for(i32 i=0; i<20; i++)
	{
		if(hasConverged(roots[0], tol) && hasConverged(roots[1], tol))
		{
			// Output roots
			// roots[0][1], roots[1][1]
			// std::cout << "AAA\n";
			zeros.push_back(roots[0][1]);
			zeros.push_back(roots[1][1]);
		
			p = p_q;
		
			return true;
		}

		divPoly(p, sigma, p_q, p_r);

		b = p_r[p_r.power()];
		a = p_r[p_r.power() - 1] - b * sigma[sigma.power() - 1];

		complex x0(0,0), x1(0,0);

		findQuadraticRoots(sigma, x0, x1);
	
		// std::cout << "quadratic roots\n";
		// std::cout << x0.real << " + " << x0.imag << "i\n";
		// std::cout << x1.real << " + " << x1.imag << "i\n";
	
		if(fabs(fabs(x0.real) - fabs(x1.real)) > 0.01 * fabs(x1.real))
		{
			// std::cout << "B\n";
			return linearShift(zeros, p, K, x, M, L, qflag, lflag, tol);
		}

		px = fabs(a - x0.real * b) + fabs(x0.imag * b);

		if(!fixed_shift && fabs(sigma[0] - tmp) / sigma[0] < step && p_x > px)
		{
			fixed_shift = true;
			fixedShiftKPoly(p, sigma, K, x0, a, b, c, d, M);
		}

		divPoly(K, sigma, k_q, k_r);

		d = k_r[k_r.power()];
		c = k_r[k_r.power() - 1] - d * sigma[sigma.power() - 1];

		tmp = sigma[sigma.power() - 2];

		sigma = nextSigma(p, sigma, K, a, b, c, d);

		K = nextKQuaraticShift(K, sigma, p_q, k_q, a, b, c, d);
		K = K.normalized();

		p_x = px;

		roots[0].push_back(x0);
		roots[1].push_back(x1);

		if(roots[0].size() > 3)
		{
			roots[0].erase(roots[0].begin());
			roots[1].erase(roots[1].begin());
		}
	}

	qflag = true;
	// std::cout << "C\n";

	return linearShift(zeros, p, K, x, M, L, qflag, lflag, tol);
}
bool jenkinsTraubPhaseThree(std::vector<complex>& roots, Polynomial& p, Polynomial& K, complex& x, u64 M, u64 L, u32 cnv, bool& quadratic_flag, bool& linear_flag, f32 tol)
{
	/*
	 * Stage 3: Find the root with variable shift iterations on the K-polynomial.
	 */

	if(cnv == 1)
	{
		// quadratic
		// std::cout << "quadratic\n";
		return quadraShift(roots, p, K, x, M, L, quadratic_flag, linear_flag, tol);
	}

	if(cnv == 2)
	{
		// linear
		// std::cout << "linear\n";
		return linearShift(roots, p, K, x, M, L, quadratic_flag, linear_flag, tol);
	}

	return false;

	// Polynomial x, r;
	// Polynomial R, E, K = H.normalized();

	// s = s - p.eval(s)/K.eval(s);
	// f32 l, s_prev;

	// for(int i=0; i<1000; i++)
	// {
	// 	if(fabs(p.eval(s)) < tol)
	// 		break;

	// 	l = -H.eval(s)/p.eval(s);
	// 	K = H + p*l;

	// 	x = Polynomial(1, {-s,1});

	// 	E = K/x;

	// 	R = E.normalized();

	// 	s_prev = s;

	// 	s = s - p.eval(s)/R.eval(s);
	// 	H = E;
	// }
}

void jenkinsTraub(std::vector<complex>& roots, Polynomial& p, f32 tol)
{
	Polynomial K;

	complex x(0.0, 0.0);

	u64 M = 20;
	u64 L = 100;

	jenkinsTraubPhaseOne(roots, p, K, x, M, tol);

	std::cout << "phase1 K:\n";
	printPoly(K);

	int cnv = jenkinsTraubPhaseTwo(roots, p, K, x, L, tol);

	if(cnv == 0)
	{
		printf("ERROR: didn't converge on phase 2\n");
		exit(1);
	}

	bool qflag = false;
	bool lflag = false;

	jenkinsTraubPhaseThree(roots, p, K, x, M, L, cnv, qflag, lflag, tol);

}

Polynomial Polynomial::root(std::vector<complex>& roots, f32 tol)
{
	Polynomial p = this->normalized();
	printf("P:\n");
	printPoly(p);
	jenkinsTraub(roots, p, tol);
	return p;
}

void Polynomial::roots(std::vector<complex>& roots, f32 tol)
{
	Polynomial r, q, p = Polynomial(*this);

	i64 n = p.power();

	for(i64 i=0; i<n; i++)
	{
		//TODO: remove leading zeros
		//TODO: normalize polynomial
		//TODO: remove zero roots

		p = p.root(roots, tol);
		// std::cout << "deflated:\n";
		// printPoly(p);
		// if(std::isnan(x))
		// 	break;
		// std::cout << "roots:\n";
		// for(complex r : roots)
		// {
		// 	std::cout << r.real << " + " << r.imag << "i\n";
		// }
		std::cout << "=================================\n";
	}
}

}
