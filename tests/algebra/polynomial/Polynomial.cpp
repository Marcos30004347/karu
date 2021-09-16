#include <assert.h>
#include <math.h>
#include "algebra/polynomial/Polynomial.hpp"

using namespace karu::algebra;

void polynomialDivisionTests()
{
	Polynomial A(3, {1,2,3,4});
	Polynomial B(2, {1,2,1});
	Polynomial Q(3,nullptr);
	Polynomial R(3,nullptr);	

	divPoly(A, B, Q, R);

	assert(Q == Polynomial(1, {-5, 4}));
	assert(R == Polynomial(1, {6, 8}));

	Q = Polynomial(0,{0});
	R = Polynomial(0,{0});
	
	Polynomial C(3, {-3,10,-5,3});
	Polynomial D(1, {1,3});

	divPoly(C, D, Q, R);

	assert(Q == Polynomial(2, {4, -2, 1}));
	assert(R == Polynomial(0, {-7}));
}

void polynomialMultiplicationTests()
{
	Polynomial A(3, {5,0,10,6});
	Polynomial B(2, {1,2,4});

	Polynomial R = mulPoly(A, B);

	assert(R == Polynomial(5, {5, 10, 30, 26, 52, 24}));
}

void polynomialAdditionTests()
{
	Polynomial A(3, {6,10,0,5});
	Polynomial B(3, {6,10,0,5});

	Polynomial C = addPoly(A,B);

	assert(C == Polynomial(3, {12, 20, 0, 10}));
}

void polynomialSubtractionTests()
{
	Polynomial A(3, {3,10,0,5});
	Polynomial B(3, {1,7,1,5});

	Polynomial C = subPoly(A,B);

	assert(C == Polynomial(2, {2, 3, -1}));
}

void polynomialRootsTests()
{
	karu::f32 tol = 0.00009;
	// std::vector<karu::f32> A_roots;
	// Polynomial A(3, {3,0,-4,1});

	// A.roots(A_roots, tol);
	// printPoly(A);

	// for(karu::f32 root : A_roots)
	// {
	// 	std::cout << root << "\n";
	// 	assert(fabs(A.eval(root)) <= tol);
	// }

	// std::vector<karu::f32> B_roots;

	// Polynomial B(2, {-3,0,2});

	// B.roots(B_roots, tol);
	// printPoly(B);
	// for(karu::f32 root : B_roots)
	// {
	// 	std::cout << root << "\n";
	// 	assert(fabs(B.eval(root)) <= tol);
	// }

	std::vector<complex> C_roots;

	// -3 + x + 3x^3 + 5x^4
	Polynomial C(4, { -3, 1, 0, 3, 5});

	std::cout << C.eval(complex(-0.07877, -0.85693)).real << "\n";
	std::cout << C.eval(complex(-0.07877, -0.85693)).real << "\n";
	std::cout << C.eval(0.70569) << "\n";

	// std::cout << (complex(1, 3) * complex(2, 1)).real << "\n";
	// std::cout << (complex(1, 3) * complex(2, 1)).imag << "\n";
	// std::cout << (complex(2, 5) * complex(4, -3)).real << "\n";
	// std::cout << (complex(2, 5) * complex(4, -3)).imag << "\n";

	printPoly(C);
	C.roots(C_roots, tol);
	printPoly(C);

	for(complex root : C_roots)
	{
		std::cout << root.real << "+" << root.imag << "i\n";
		assert(fabs(C.eval(root)) <= tol);
	}
}

int main()
{
	polynomialDivisionTests();
	polynomialMultiplicationTests();
	polynomialAdditionTests();
	polynomialSubtractionTests();
	polynomialRootsTests();

	return 0;
}
