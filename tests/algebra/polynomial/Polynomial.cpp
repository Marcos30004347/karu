#include <assert.h>
#include <math.h>
#include "algebra/polynomial/Polynomial.hpp"

using namespace karu::algebra;

void polynomialDivisionTests()
{
	RealPoly A(3, {1,2,3,4});
	RealPoly B(2, {1,2,1});
	RealPoly Q(3,nullptr);
	RealPoly R(3,nullptr);	

	divPoly(A, B, Q, R);

	assert(Q == RealPoly(1, {-5, 4}));
	assert(R == RealPoly(1, {6, 8}));

	Q = RealPoly(0,{0});
	R = RealPoly(0,{0});
	
	RealPoly C(3, {-3,10,-5,3});
	RealPoly D(1, {1,3});

	divPoly(C, D, Q, R);

	assert(Q == RealPoly(2, {4, -2, 1}));
	assert(R == RealPoly(0, {-7}));
}

void polynomialMultiplicationTests()
{
	RealPoly A(3, {5,0,10,6});
	RealPoly B(2, {1,2,4});

	RealPoly R = mulPoly(A, B);

	assert(R == RealPoly(5, {5, 10, 30, 26, 52, 24}));
}

void polynomialAdditionTests()
{
	RealPoly A(3, {6,10,0,5});
	RealPoly B(3, {6,10,0,5});

	RealPoly C = addPoly(A,B);

	assert(C == RealPoly(3, {12, 20, 0, 10}));
}

void polynomialSubtractionTests()
{
	RealPoly A(3, {3,10,0,5});
	RealPoly B(3, {1,7,1,5});

	RealPoly C = subPoly(A,B);

	assert(C == RealPoly(2, {2, 3, -1}));
}

void polynomialRootsTests()
{
	// std::vector<karu::f32> A_roots;
	// RealPoly A(3, {3,0,-4,1});

	// A.roots(A_roots, tol);
	// printPoly(A);

	// for(karu::f32 root : A_roots)
	// {
	// 	std::cout << root << "\n";
	// 	assert(fabs(A.eval(root)) <= tol);
	// }

	// std::vector<karu::f32> B_roots;

	// RealPoly B(2, {-3,0,2});

	// B.roots(B_roots, tol);
	// printPoly(B);
	// for(karu::f32 root : B_roots)
	// {
	// 	std::cout << root << "\n";
	// 	assert(fabs(B.eval(root)) <= tol);
	// }


	// -3 + x + 3x^3 + 5x^4

	// std::cout << C.eval(complex(-0.07877, -0.85693)).real << "\n";
	// std::cout << C.eval(complex(-0.07877, -0.85693)).real << "\n";
	// std::cout << C.eval(0.70569) << "\n";

	// std::cout << (complex(1, 3) * complex(2, 1)).real << "\n";
	// std::cout << (complex(1, 3) * complex(2, 1)).imag << "\n";
	// std::cout << (complex(2, 5) * complex(4, -3)).real << "\n";
	// std::cout << (complex(2, 5) * complex(4, -3)).imag << "\n";

	RealPoly C(4, { -3, 1, 0, 3, 5});

	std::vector<complex> C_roots = polyRoots(C);

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
