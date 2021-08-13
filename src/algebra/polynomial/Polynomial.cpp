#include "algebra/polynomial/Polynomial.hpp"
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>

namespace karu::algebra {

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

f32 Polynomial::eval(f32 x)
{
	f32 val = 0;

	for(i64 i=0; i<=this->power(); i++)
	{
		val += this->p_coefficients[i]*pow(x, i);
	}

	return val;
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

f32 newtonRaphson(Polynomial& p, f32 tol)
{
	Polynomial pdx = p.dx();

	f32 x = 0.0;

	while(pdx.eval(x) == 0)
	{
		x += randomBetweenZeroOne();
	}

	while(fabs(p.eval(x)) > tol)
	{
		x = x - p.eval(x)/pdx.eval(x);
	}

	return x;
}


void jenkinsTraubPhaseOne(Polynomial& p, Polynomial& H, f32& s, u64 M, f32 tol)
{
	f32 l;
	
	Polynomial K, r, x = Polynomial(1, {0,1});
	
	H = p.dx();
	
	s = 0.0;
	
	for(i64 i=0; i<M; i++)
	{
		l = -H.eval(0)/p.eval(0);
		K = H + p*l;
		H = K/x;
	}
}

void jenkinsTraubPhaseTwo(Polynomial& p, Polynomial& H, f32& s, u64 L, f32 tol)
{
	Polynomial K, E, U, O, r, x;
	Polynomial cauchy = p.cauchy();
	f32 l, t_curr, t_prev, t_next;
	bool stop = false;

	u64 it = 0;
	do
	{
		s = newtonRaphson(cauchy, tol);
	
		if(p.eval(s) == 0)
			break;
	
		for(i64 i=0; i<L && p.eval(s) != 0; i++)
		{
			l = -H.eval(s)/p.eval(s);
			K = H + p*l;

			x = Polynomial(1, {-s,1});

			E = K/x;

			U = H.normalized();
			O = E.normalized();

			t_curr = s - p.eval(s)/U.eval(s);
			t_next = s - p.eval(s)/O.eval(s);

			if(
				i > 0 &&
				fabs(t_curr - t_prev)<= 0.5*fabs(t_prev) &&
				fabs(t_next - t_curr)<=0.5*fabs(t_curr)
			) {
				stop = true;
				break;
			}
			
			t_prev = t_curr;

			H = E;
		}
	} while(!stop && it++ < 1000);
}

void jenkinsTraubPhaseThree(Polynomial& p, Polynomial& H, f32& s, u64 L, f32 tol)
{
	Polynomial x, r;
	Polynomial R, E, K = H.normalized();
	
	s = s - p.eval(s)/K.eval(s);
	f32 l, s_prev;

	for(int i=0; i<1000; i++)
	{
		if(fabs(p.eval(s)) < tol)
			break;
		
		l = -H.eval(s)/p.eval(s);
		K = H + p*l;

		x = Polynomial(1, {-s,1});

		E = K/x;

		R = E.normalized();
	
		s_prev = s;
	
		s = s - p.eval(s)/R.eval(s);
		H = E;
	}
}

f32 jenkinsTraub(Polynomial p, f32 tol)
{
	Polynomial H;
	
	f32 s = 0.0;
	u64 M = 5, L = 100;

	jenkinsTraubPhaseOne(p, H, s, M, tol);
	jenkinsTraubPhaseTwo(p, H, s, L, tol);
	jenkinsTraubPhaseThree(p, H, s, L, tol);

	return s;
}

f32 Polynomial::root(f32 tol)
{
	if(this->power() == 1)
		return -1*this->p_coefficients[0]/this->p_coefficients[1];

	return jenkinsTraub(this->normalized(), tol);
}

void Polynomial::roots(std::vector<f32>& roots, f32 tol)
{
	Polynomial r, q, p = Polynomial(*this);

	f32 x;

	i64 n = p.power();

	for(i64 i=0; i<n; i++)
	{
		if(p.power() == 1)
		{
			x = -1*p.p_coefficients[0]/p.p_coefficients[1];
			
			if(!std::isnan(x))
				roots.push_back(x);
			
			break;
		}
		if(p.power() == 2)
		{
			f32 c = p[0], b = p[1], a = p[2];
	
			f32 delta = b*b - 4*a*c;

			x = (-b+sqrt(delta))/(2*a);
	
			if(!std::isnan(x))
				roots.push_back(x);		
			
			x = (-b-sqrt(delta))/(2*a);

			if(!std::isnan(x))
				roots.push_back(x);		

			break;
		}
		else
		{
			x = p.root(tol);

			if(std::isnan(x))
				break;

			roots.push_back(x);

			q = Polynomial(1, {-x,1});
			p = p/q;
		}
	}
}

}
