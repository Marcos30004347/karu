#include "algebra/polynomial/Polynomial.hpp"

#include <cmath>
#include <random>
#include <chrono>

namespace karu::algebra {

complex::complex(f32 r, f32 i)
{
	this->real = r;
	this->imag = i;
}

complex::complex(f32 r)
{
	this->real = r;
	this->imag = 0;
}

complex::complex()
{
	this->real = 0;
	this->imag = 0;
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

	complex n = a * other.conj();
	complex d = b * other.conj();

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
	return complex(this->real + other, this->imag);
}

complex complex::operator-(f32 other)
{
	return complex(this->real - other, this->imag);
}

complex complex::operator/(f32 other)
{
	return complex(this->real / other, this->imag / other);
}

void printPoly(RealPoly& p)
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

void divPoly(const RealPoly p, const RealPoly& d, RealPoly& q, RealPoly& r)
{
	i64 i,j;
	i64 power = p.power() - d.power();
	f32 ratio;

	if(power < 0) return;

	q = RealPoly(power, {});

	r = RealPoly(p);

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


void divPoly(const RealPoly p, const f32 d, RealPoly& q)
{
	q = RealPoly(p);
	for (int j=0; j<=p.power(); j++)
		q[j] = p[j]/d;
}

RealPoly addPoly(const RealPoly& A, const RealPoly& B)
{
	RealPoly X = A.power() > B.power() ? RealPoly(A) : RealPoly(B);
	RealPoly Y = A.power() > B.power() ? RealPoly(B) : RealPoly(A);

	for(i64 i=0; i<=Y.power(); i++)
	{
		X[i] += Y[i];
	}

	while(X[X.p_power] == 0) { X.p_power--; }
	return X;
}

RealPoly subPoly(const RealPoly& A, const RealPoly& B)
{
	RealPoly B_ = mulPoly(B, -1.0);

	RealPoly X = A.power() > B.power() ? RealPoly(A) : RealPoly(B_);
	RealPoly Y = A.power() > B.power() ? RealPoly(B_) : RealPoly(A);

	for(i64 i=0; i<=Y.power(); i++)
	{
		X[i] += Y[i];
	}

	while(X[X.p_power] == 0) { X.p_power--; }

	return X;
}

RealPoly mulPoly(const RealPoly& A, const RealPoly& B)
{
	u64 m = A.power() + 1;
	u64 n = B.power() + 1;

	RealPoly AB(m+n-2,{});

	for (int i=0; i<m; i++)
	{
		for (int j=0; j<n; j++)
		{
			AB[i+j] += A[i]*B[j];
		}
	}
	return AB;
}

RealPoly mulPoly(const RealPoly& A, f32 B)
{
	RealPoly AB(A);
	for (int j=0; j<=AB.power(); j++)
		AB[j] = A[j]*B;
	return AB;
}

RealPoly::RealPoly()
{
	this->p_power = 0;
	this->p_coefficients = new f32[this->p_power+1];
	this->p_coefficients[0] = 0.0;
}

RealPoly::RealPoly(const RealPoly& other)
{
	this->p_power = other.p_power;
	this->p_coefficients = new f32[this->p_power+1];

	std::copy(other.p_coefficients, other.p_coefficients + other.p_power+1, this->p_coefficients);
}

RealPoly::RealPoly(u64 power, std::initializer_list<f32> coefs)
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

RealPoly::RealPoly(u64 power)
{
	i64 i;
	this->p_power = power;
	this->p_coefficients = new f32[++power];
	std::fill(&this->p_coefficients[0], &this->p_coefficients[power-1]+1, 0);
}


RealPoly::RealPoly(u64 power, f32* coefs)
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

RealPoly::~RealPoly(){
	delete this->p_coefficients;
}

i64 RealPoly::power() const
{
	return this->p_power;
}


RealPoly RealPoly::operator/(const RealPoly& other)
{
	RealPoly q,r;
	divPoly(*this, other, q, r);
	return q;
}

RealPoly RealPoly::operator/(const f32 other)
{
	RealPoly q;
	divPoly(*this, other, q);
	return q;
}

RealPoly RealPoly::operator*(const f32 other)
{
	return mulPoly(*this, other);
}

RealPoly RealPoly::operator*(const RealPoly& other)
{
	return mulPoly(*this, other);
}

RealPoly RealPoly::operator+(const RealPoly& other)
{
	return addPoly(*this, other);
}

RealPoly RealPoly::operator-(const RealPoly& other)
{
	return subPoly(*this, other);
}

RealPoly& RealPoly::operator = (const RealPoly& other)
{
	this->p_power = other.p_power;
	this->p_coefficients = new f32[this->p_power+1];
	std::copy(other.p_coefficients, other.p_coefficients + other.p_power+1, this->p_coefficients);
	return *this;
}

bool RealPoly::operator == (const RealPoly& other)
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

f32& RealPoly::operator[](u64 i) const
{
	return this->p_coefficients[i];
}

// horners method
f32 RealPoly::eval(f32 x)
{
	f32 px = this->p_coefficients[this->power()];

	for(i64 i = this->power() - 1; i >= 0; i--)
	{
		px = px * x + this->p_coefficients[i];
	}

	return px;
}

f32 syntheticDivision(RealPoly& p, f32 x, RealPoly& q)
{
	q = RealPoly(p.power() - 1, {});
	
	q[p.power() - 1] = p[p.power()];

	for(i32 i = p.power() - 2; i >= 0; i--)
	{
		q[i] = p[i + 1] + q[i + 1] * x;
	}

	return p[0] + q[0] * x;
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
complex RealPoly::eval(complex x)
{
	complex px(0, 0);

	for(i64 i = this->power(); i >= 0; i--)
	{
		px = px * x + this->p_coefficients[i];
	}

	return px;
}

RealPoly RealPoly::dx()
{
	RealPoly d(this->power() - 1, {});

	for(i64 i=1; i<=this->power(); i++)
	{
		d[i-1] = this->p_coefficients[i]*i;
	}

	return d;
}

RealPoly RealPoly::abs()
{
	RealPoly d(this->power(), {});

	for(i64 i=0; i<=this->power(); i++)
	{
		d[i] = this->p_coefficients[i] >= 0 ? this->p_coefficients[i] : -this->p_coefficients[i];
	}

	return d;
}

RealPoly RealPoly::reverse()
{
	RealPoly d(this->power(), {});

	for(i64 i=0; i<=this->power(); i++)
	{
		d[this->power() - i] = this->p_coefficients[i];
	}

	return d;
}

RealPoly RealPoly::noLeadingTerm()
{
	RealPoly d(this->power() - 1, {});

	for(i64 i=0; i<this->power(); i++)
	{
		d[i] = this->p_coefficients[i];
	}


	return d;
}

RealPoly RealPoly::noConstantTerm()
{
	RealPoly d(this->power(), {});

	for(i64 i=0; i<=this->power(); i++)
	{
		d[i] = this->p_coefficients[i];
	}

	d[0] = 0;

	return d;
}

RealPoly RealPoly::normalized()
{

	RealPoly n(this->power(), {});

	if(n[n.power()] == 1)
		return n;

	for(i64 i=0; i<=this->power(); i++)
	{
		n[i] = this->p_coefficients[i]/this->p_coefficients[this->p_power];
	}

	return n;
}

RealPoly RealPoly::cauchy()
{
	RealPoly p = this->normalized();

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

f32 newtonRaphson(RealPoly& p, f32 x0, f32 tol, u32 max_it = 300)
{
	RealPoly pdx = p.dx();

	f32 x = x0;
	f32 prev = 0;
	for(i32 i = 0; i<max_it; i++)
	{
		prev = x;
	
		x = x - p.eval(x)/pdx.eval(x);

    if (fabs(prev - x) < std::numeric_limits<f32>::epsilon()) 
		{
			break;
		}
	}

	return x;
}

f32 rootRadius(RealPoly& p, f32 tol = 0.01)
{
	RealPoly p_ = p.abs();

	p_[0] *= -1;

	printPoly(p_);

	return newtonRaphson(p_, 1.0, 0.01, 100);
}

void jenkinsTraubPhaseOne(std::vector<complex>& roots, RealPoly& p, RealPoly& K, complex& s, u64 M, f32 tol)
{
	RealPoly z = RealPoly(1, { 0,1 });

	K = p.dx() / (p.power() + 1);

	for(i64 i=1; i<M; i++)
	{
		K = (K.noConstantTerm() + p.noConstantTerm()*(-K[0]/p[0]))/z;
	}
}

RealPoly nextSigma(RealPoly& p, RealPoly& sigma, RealPoly& K, f32& a, f32& b, f32& c, f32& d)
{
	f32 u = sigma[1];
	f32 v = sigma[0];

	f32 b1 = -K[0] / p[0];
	f32 b2 = -(K[1] + b1 * p[1]) / p[0];

	f32 a1 = b * c - a * d;
	f32 a2 = a * c + u * a * d + v * b * d;
	f32 c2 = b1 * a2;
	f32 c3 = b1 * b1 * (a * a + u * a * b + v * b * b);
	f32 c4 = v * b2 * a1 - c2 - c3;
	f32 c1 = c * c + u * c * d + v * d * d + b1 * (a * c + u * b * c + v * b * d) - c4;
	f32 du = -(u * (c2 + c3) + v * (b1 * a1 + b2 * a2)) / c1;
	f32 dv = v * c4 / c1;

	return RealPoly(2, {v + dv, u + du, 1});
}

RealPoly nextKQuaraticShift(RealPoly K, RealPoly& sigma, RealPoly& p_q, RealPoly& k_q, f32& a, f32& b, f32& c, f32& d)
{
	f32 t = (a * a + sigma[1] * a * b + sigma[0] * b * b) / (b * c - a * d);
  
	printPoly(sigma);

	RealPoly lin(1, {});

  lin[0] = -(a * c + sigma[1] * a * d + sigma[0] * b * d) / (b * c - a * d);
	lin[1] = 1;
	
	K = (k_q * t) + (lin * p_q);

	K[0] += b;

	return K;
}

int fixedShiftKPoly(RealPoly& p, RealPoly& sigma, RealPoly& K, complex& x, f32& a, f32& b, f32& c, f32& d, u32 max_it)
{
	RealPoly p_r, p_q, k_q, k_r;

	sigma[0] = x.real*x.real + x.imag*x.imag;
	sigma[1] = -2.0 * x.real;
	sigma[2] = 1;

	divPoly(p, sigma, p_q, p_r);

	b = p_r[p_r.power()];
	a = p_r[1] - b * sigma[1];

	complex p_x = x.conj() * -b + a;

	std::cout << p_x.real << " + " << p_x.imag << "i\n";

	complex t_l[3] = {complex(0,0),complex(0,0),complex(0,0)};
	f32 s_l[3] 		 = {0,0,0};

	for(i32 j=0; j < max_it; j++)
	{
    std::cout << "*********\n";

		printPoly(K);
		K = K.normalized();
		printPoly(K);

		divPoly(K, sigma, k_q, k_r);

		d = k_r[k_r.power()];
		c = k_r[k_r.power() - 1] - d * sigma[sigma.power() - 1];
	
    std::cout << d << " " << c << "\n";

		RealPoly s = nextSigma(p, sigma, K, a, b, c, d);

		complex k_x = x.conj() * -d + c;

		t_l[0] = t_l[1];
		t_l[1] = t_l[2];
	
		s_l[0] = s_l[1];
		s_l[1] = s_l[2];

		t_l[2] = x - p_x / k_x;
		s_l[2] = s[0];
	
		std::cout << t_l[0].real << " + " << t_l[0].imag << "i, ";
		std::cout << t_l[1].real << " + " << t_l[1].imag << "i, ";
		std::cout << t_l[2].real << " + " << t_l[2].imag << "i\n";
	
		std::cout << s_l[0] << " " << s_l[1] << " " << s_l[2] << "\n";

		if(
			std::fabs(s_l[1] - s_l[0]) < std::fabs(s_l[0]) / 2.0 &&
			std::fabs(s_l[2] - s_l[1]) < std::fabs(s_l[1]) / 2.0
		) return 2;

		if(
			(t_l[1] - t_l[0]).abs() < (t_l[0]).abs() / 2.0 &&
			(t_l[2] - t_l[1]).abs() < (t_l[1]).abs() / 2.0
		) return 1;

		K = nextKQuaraticShift(K, sigma, p_q, k_q, a, b, c, d);
	}

	return 0;
}


int jenkinsTraubPhaseTwo(std::vector<complex>& roots, RealPoly& p, RealPoly& K, complex& x, u64 L, f32 tol)
{
	u32 i, j;

	f32 a, b, d, c;

	f32 deg2rad = M_PI / 180.0;
  f32 phi = 49.0 * deg2rad;

	f32 radius = rootRadius(p, 1e-2);

	std::cout << "root_radius" << "\n";
	std::cout << radius << "\n";

	RealPoly sigma(2, {});

	for(i=0; i < 20; i++)
	{
		x = complex(radius, 0) * complex(cos(phi), sin(phi));
	
		std::cout << "root" << "\n";
		std::cout << x.real << " + " << x.imag << "i\n";
	
		i32 j = fixedShiftKPoly(p, sigma, K, x, a, b, c, d, 20 * (i + 1));
		
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
void findQuadraticRoots(RealPoly& p, complex& x1, complex& x2)
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


bool linearShift(std::vector<complex>& zeros, RealPoly& p, RealPoly& K, complex& x, u64 M, u64 L, bool& qflag, bool& lflag, f32 tol)
{
	if(lflag)
	{
		return false;
	}

	std::vector<f32> roots;

	RealPoly defl_p, defl_k;
	
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

		px = syntheticDivision(p, s, defl_p);
		// defl_p = p / RealPoly(1, {-s, 1});
		// px = p.eval(s);

		if(fabs(px) <= tol)
		{
			// std::cout << "2#\n";
	
			zeros.push_back(complex(s, 0));
			p = defl_p;
	
			return true;
		}
	
		// defl_k = K / RealPoly(1, { -s, 1 });
		kx = syntheticDivision(K, s, defl_k);
		// kx = K.eval(s);
		
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

bool quadraShift(std::vector<complex>& zeros, RealPoly& p, RealPoly& K, complex& x, u64& M, u64& L, bool& qflag, bool& lflag, f32 tol)
{
	if(qflag)
	{
		return false;
	}

	complex x0, x1;

	f32 step = 0.01;
	f32 a, b, c, d;
	f32 px = 0, p_x = 0, tmp = 0;

	std::vector<complex> roots[2];
	
	RealPoly p_q, p_r, k_q, k_r, sigma(2, {});

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
			zeros.push_back(roots[0][1]);
			zeros.push_back(roots[1][1]);
		
			p = p_q;
		
			return true;
		}

		divPoly(p, sigma, p_q, p_r);

		b = p_r[p_r.power()];
		a = p_r[p_r.power() - 1] - b * sigma[sigma.power() - 1];

	
    std::cout << "quadratic shift sigma:\n";
		printPoly(sigma);

		findQuadraticRoots(sigma, x0, x1);

		if(fabs(fabs(x0.real) - fabs(x1.real)) > 0.01 * fabs(x1.real))
		{
			std::cout << "applying linearShift\n";
			return linearShift(zeros, p, K, x, M, L, qflag, lflag, tol);
		}

		px = fabs(a - x0.real * b) + fabs(x0.imag * b);
	
		std::cout << "px: " << px << "\n";
	
		if(!fixed_shift && fabs(sigma[0] - tmp) / sigma[0] < step && p_x > px)
		{
			fixed_shift = true;
      std::cout << "apply fixed shift\n";
			fixedShiftKPoly(p, sigma, K, x0, a, b, c, d, M);
		}

		divPoly(K, sigma, k_q, k_r);

		d = k_r[k_r.power()];
		c = k_r[k_r.power() - 1] - d * sigma[sigma.power() - 1];
	
		std::cout << "d" << " " << "c" << "\n";
		std::cout << d << " " << c << "\n";
	
		tmp = sigma[0];
		sigma = nextSigma(p, sigma, K, a, b, c, d);
		std::cout << "sigma\n";
		printPoly(sigma);

		K = nextKQuaraticShift(K, sigma, p_q, k_q, a, b, c, d);
		K = K.normalized();
		std::cout << "K\n";
		printPoly(K);
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
bool jenkinsTraubPhaseThree(std::vector<complex>& roots, RealPoly& p, RealPoly& K, complex& x, u64 M, u64 L, u32 cnv, bool& quadratic_flag, bool& linear_flag, f32 tol)
{
	/*
	 * Stage 3: Find the root with variable shift iterations on the K-polynomial.
	 */
  
	std::cout << "phase3 K:\n";
	printPoly(K);
	std::cout << cnv << "\n";

	if(cnv == 2)
	{
		// quadratic
		// std::cout << "quadratic\n";
		return quadraShift(roots, p, K, x, M, L, quadratic_flag, linear_flag, tol);
	}

	if(cnv == 1)
	{
		// linear
		// std::cout << "linear\n";
		return linearShift(roots, p, K, x, M, L, quadratic_flag, linear_flag, tol);
	}

	return false;

	// RealPoly x, r;
	// RealPoly R, E, K = H.normalized();

	// s = s - p.eval(s)/K.eval(s);
	// f32 l, s_prev;

	// for(int i=0; i<1000; i++)
	// {
	// 	if(fabs(p.eval(s)) < tol)
	// 		break;

	// 	l = -H.eval(s)/p.eval(s);
	// 	K = H + p*l;

	// 	x = RealPoly(1, {-s,1});

	// 	E = K/x;

	// 	R = E.normalized();

	// 	s_prev = s;

	// 	s = s - p.eval(s)/R.eval(s);
	// 	H = E;
	// }
}

void findLinearRoots(RealPoly& p, complex& x)
{
	x = complex(-p[0]/p[1], 0);
}

void jenkinsTraub(std::vector<complex>& roots, RealPoly& p, f32 tol)
{
	if(p.power() == 0)
	{
		return;
	}

	if(p.power() == 1)
	{
		complex x;
		findLinearRoots(p, x);
		roots.push_back(x);
		return;
	}

	if(p.power() == 2)
	{
		complex x1, x2;
		findQuadraticRoots(p, x1, x2);
		roots.push_back(x1);
		roots.push_back(x2);
		return;
	}

	RealPoly K;

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

RealPoly RealPoly::root(std::vector<complex>& roots, f32 tol)
{
	RealPoly p = this->normalized();
	jenkinsTraub(roots, p, tol);
	return p;
}

RealPoly removeZeroRoots(std::vector<complex>& roots, RealPoly& p)
{
	i32 i = 0;
	
	while(p[i] == 0)
	{
		i++;
		roots.push_back(complex(0,0));
	}

	RealPoly p_ = RealPoly(p.power() -  i, {});

	for(i32 j = 0; j <= p.power() - i; j++)
	{
		p_[p.power() - j] = p[p.power() - j];
	}

	return p_;
}

std::vector<complex> polyRoots(RealPoly p, f32 tol)
{
	std::vector<complex> roots;

	i64 n = p.power();

	while(roots.size() < n)
	{
		p = p.normalized();
		p = removeZeroRoots(roots, p);
		p = p.root(roots, tol);
	}

	return roots;
}

}
