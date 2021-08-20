#include <assert.h>

#include "algebra/SVD/Householder.hpp"
#include "algebra/SVD/GolubReinschSVD.hpp"

using namespace karu;
using namespace karu::algebra;

void householderTests()
{
	Matrix a(4,3, {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12,
	});

	Matrix u, v;

	f32 e[3];
	f32 q[3];

	householderBidiagonalForm(a, q, e, u, v, 2.22e-16);

	std::cout << "e:\n";
	for(int i=0; i<3; i++)
		std::cout << e[i] << ", ";
	
	std::cout << "\n";
	std::cout << "q:\n";
	
	for(int i=0; i<3; i++)
		std::cout << q[i] << ", ";
	
	std::cout << "\n";
	std::cout << "u:\n";
	printMatrix(u);
	std::cout << "v:\n";
	printMatrix(v);
	std::cout << "a:\n";

	Matrix b(4, 3, {
		q[0], e[1], 0,
		0, q[1], e[2],
		0,   0,  q[2],
		0, 0, 0,
	});

	std::cout << "u*b*v':\n";
	printMatrix(u*b*transpose(v));
	std::cout << "diff':\n";
	printMatrix((u*b*transpose(v)) - a);

}


void bidiagonalTests()
{
	Matrix a(4,3, {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12,
	});
	Matrix u, v;

	barlowBidiagonalization(a, u, v);
}


void svdTests()
{
	Matrix a(4,3, {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12,
	});
	Matrix u, v;
	f32 s[3];

	f32 eps = 1.e-15;
	f32 tol = 1.e-64/eps;

	golubReinschSVD(a, u, s, v, tol, eps);
	std::cout << "U\n";
	printMatrix(u);
	std::cout << "V'\n";
	printMatrix(v);
	for(int i=0; i<3; i++)
		std::cout << s[i] << " ";
	std::cout << "\n";
	std::cout << "\n";

	Matrix b(4, 3, {
		s[0], 0, 0,
		0, s[1], 0,
		0, 0,  s[2],
		0, 0,    0,
	});

	printMatrix(u*b*transpose(v));
	std::cout << "\n";
	printMatrix(u*b*transpose(v) - a);
	std::cout << "\n";

	f32 n = v[0][2]*v[0][2] + v[1][2]*v[1][2] + v[2][2]*v[2][2];

	n = sqrt(n);


}

int main()
{
	
	householderTests();
	svdTests();
	bidiagonalTests();
	return 0;
}
