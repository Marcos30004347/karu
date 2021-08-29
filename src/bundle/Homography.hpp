#pragma once

#include <cmath>
#include "bundle/codegen/Homography.hpp"
#include "bundle/codegen/Estimation.hpp"
#include "algebra/matrix/Matrix.hpp"
#include "algebra/polynomial/Polynomial.hpp"
#include "algebra/linear/Linear.hpp"
// #include "algebra/linear/SingularValueDecomposition.hpp"
// #include "algebra/linear/svd.hpp"
#include "algebra/SVD/SVD.hpp"
// #include "algebra/linear/Eigen/SVD"

using namespace karu;
using namespace karu::algebra;


f32 betaP(f32* b, f32* c, f32* d, f32* p11, f32* p21, f32* p12, f32* p22, f32* a)
{
	f32 beta0 = -cos(a[1])*c[0]*d[6]*p22[1]*sin(a[0]);
	// std::cout << "b: "<< a[1] << "\n";
	// std::cout << "b: "<< c[0] << "\n";
	// std::cout << "b: "<< d[6] << "\n";
	// std::cout << "b: "<< p22[1] << "\n";
	beta0 += cos(a[1])*c[3]*d[6]*p21[0]*sin(a[0]);
	beta0 += cos(a[1])*c[6]*d[0]*p22[1]*sin(a[0]);

	beta0 -= cos(a[1])*c[6]*d[3]*p21[0]*sin(a[0]);
	beta0 -= cos(a[1])*cos(a[0])*c[3]*d[6]*p21[1];
	beta0 += cos(a[1])*cos(a[0])*c[3]*d[6]*p22[1];

	beta0 += cos(a[1])*cos(a[0])*c[6]*d[3]*p21[1];
	beta0 -= cos(a[1])*cos(a[0])*c[6]*d[3]*p22[1];
	beta0 -= c[0]*d[6]*p21[0]*sin(a[1])*sin(a[0]);

	beta0 += c[0]*d[6]*p22[0]*sin(a[1])*sin(a[0]);
	beta0 += c[6]*d[0]*p21[0]*sin(a[1])*sin(a[0]);
	beta0 -= c[6]*d[0]*p22[0]*sin(a[1])*sin(a[0]);

	beta0 += cos(a[0])*c[0]*d[6]*p21[1]*sin(a[1]);
	beta0 -= cos(a[0])*c[3]*d[6]*p22[0]*sin(a[1]);
	beta0 -= cos(a[0])*c[6]*d[0]*p21[1]*sin(a[1]);

	beta0 += cos(a[0])*c[6]*d[3]*p22[0]*sin(a[1]);
	beta0 += cos(a[1])*c[0]*d[3]*sin(a[0]);
	beta0 -= cos(a[1])*c[3]*d[0]*sin(a[0]);

	beta0 -= cos(a[0])*c[0]*d[3]*sin(a[1]);
	beta0 += cos(a[0])*c[3]*d[0]*sin(a[1]);

	f32 beta1 = +cos(a[1])*b[0]*c[6]*p22[1]*sin(a[0]);
	beta1 += cos(a[1])*b[3]*c[6]*p21[0]*sin(a[0]);
	beta1 += cos(a[1])*b[6]*c[0]*p22[1]*sin(a[0]);

	beta1 -= cos(a[1])*b[6]*c[3]*p21[0]*sin(a[0]);
	beta1 -= cos(a[1])*cos(a[0])*b[3]*c[6]*p21[1];
	beta1 += cos(a[1])*cos(a[0])*b[3]*c[6]*p22[1];

	beta1 += cos(a[1])*cos(a[0])*b[6]*c[3]*p21[1];
	beta1 -= cos(a[1])*cos(a[0])*b[6]*c[3]*p22[1];
	beta1 -= b[0]*c[6]*p21[0]*sin(a[0])*sin(a[1]);

	beta1 += b[0]*c[6]*p22[0]*sin(a[0])*sin(a[1]);
	beta1 += b[6]*c[0]*p21[0]*sin(a[0])*sin(a[1]);
	beta1 -= b[6]*c[0]*p22[0]*sin(a[0])*sin(a[1]);

	beta1 += cos(a[0])*b[0]*c[6]*p21[1]*sin(a[1]);
	beta1 -= cos(a[0])*b[3]*c[6]*p22[0]*sin(a[1]);
	beta1 -= cos(a[0])*b[6]*c[0]*p21[1]*sin(a[1]);

	beta1 += cos(a[0])*b[6]*c[3]*p22[0]*sin(a[1]);
	beta1 += cos(a[1])*b[0]*c[3]*sin(a[0]);
	beta1 -= cos(a[1])*b[3]*c[0]*sin(a[0]);

	beta1 -= cos(a[0])*b[0]*c[3]*sin(a[1]);
	beta1 += cos(a[0])*b[3]*c[0]*sin(a[1]);

	return beta0/beta1;

}

f32 gammaP(f32* b, f32* c, f32* d, f32* p11, f32* p21, f32* p12, f32* p22, f32* a)
{

	f32 gama0 = -cos((a[1]))*(b[0])*(d[6])*(p22[1])*sin((a[0]));
	gama0 += cos((a[1]))*(b[3])*(d[6])*(p21[0])*sin((a[0]));
	gama0 += cos((a[1]))*(b[6])*(d[0])*(p22[1])*sin((a[0]));

	gama0 -= cos((a[1]))*(b[6])*(d[3])*(p21[0])*sin((a[0]));
	gama0 -= cos((a[1]))*cos((a[0]))*(b[3])*(d[6])*(p21[1]);
	gama0 += cos((a[1]))*cos((a[0]))*(b[3])*(d[6])*(p22[1]);

	gama0 += cos((a[1]))*cos((a[0]))*(b[6])*(d[3])*(p21[1]);
	gama0 -= cos((a[1]))*cos((a[0]))*(b[6])*(d[3])*(p22[1]);
	gama0 -= (b[0])*(d[6])*(p21[0])*sin((a[0]))*sin((a[1]));

	gama0 += (b[0])*(d[6])*(p22[0])*sin((a[0]))*sin((a[1]));
	gama0 += (b[6])*(d[0])*(p21[0])*sin((a[0]))*sin((a[1]));
	gama0 -= (b[6])*(d[0])*(p22[0])*sin((a[0]))*sin((a[1]));

	gama0 += cos((a[0]))*(b[0])*(d[6])*(p21[1])*sin((a[1]));
	gama0 -= cos((a[0]))*(b[3])*(d[6])*(p22[0])*sin((a[1]));
	gama0 -= cos((a[0]))*(b[6])*(d[0])*(p21[1])*sin((a[1]));

	gama0 += cos((a[0]))*(b[6])*(d[3])*(p22[0])*sin((a[1]));
	gama0 += cos((a[1]))*(b[0])*(d[3])*sin((a[0]));
	gama0 -= cos((a[1]))*(b[3])*(d[0])*sin((a[0]));

	gama0 -= cos((a[0]))*(b[0])*(d[3])*sin((a[1]));
	gama0 += cos((a[0]))*(b[3])*(d[0])*sin((a[1]));

	f32 gama1 = -cos((a[1]))*(b[0])*(c[6])*(p22[1])*sin((a[0]));
	gama1 += cos((a[1]))*(b[3])*(c[6])*(p21[0])*sin((a[0]));
	gama1 += cos((a[1]))*(b[6])*(c[0])*(p22[1])*sin((a[0]));

	gama1 -= cos((a[1]))*(b[6])*(c[3])*(p21[0])*sin((a[0]));
	gama1 -= cos((a[1]))*cos((a[0]))*(b[3])*(c[6])*(p21[1]);
	gama1 += cos((a[1]))*cos((a[0]))*(b[3])*(c[6])*(p22[1]);

	gama1 += cos((a[1]))*cos((a[0]))*(b[6])*(c[3])*(p21[1]);
	gama1 -= cos((a[1]))*cos((a[0]))*(b[6])*(c[3])*(p22[1]);
	gama1 -= (b[0])*(c[6])*(p21[0])*sin((a[0]))*sin((a[1]));

	gama1 += (b[0])*(c[6])*(p22[0])*sin((a[0]))*sin((a[1]));
	gama1 += (b[6])*(c[0])*(p21[0])*sin((a[0]))*sin((a[1]));
	gama1 -= (b[6])*(c[0])*(p22[0])*sin((a[0]))*sin((a[1]));

	gama1 += cos((a[0]))*(b[0])*(c[6])*(p21[1])*sin((a[1]));
	gama1 -= cos((a[0]))*(b[3])*(c[6])*(p22[0])*sin((a[1]));
	gama1 -= cos((a[0]))*(b[6])*(c[0])*(p21[1])*sin((a[1]));

	gama1 += cos((a[0]))*(b[6])*(c[3])*(p22[0])*sin((a[1]));
	gama1 += cos((a[1]))*(b[0])*(c[3])*sin((a[0]));
	gama1 -= cos((a[1]))*(b[3])*(c[0])*sin((a[0]));

	gama1 -= cos((a[0]))*(b[0])*(c[3])*sin((a[1]));
	gama1 += cos((a[0]))*(b[3])*(c[0])*sin((a[1]));

	return -gama0 / gama1;
}

// Matrix coplanarThreePointHomography(Matrix image0[3], Matrix image1[3], f32 angles[3])
// {
// 	f32 p1[2] = { image0[0][0][0], image0[0][1][0] };
// 	f32 p2[2] = { image1[0][0][0], image1[0][1][0] };

// 	f32 p3[2] = { image0[1][0][0], image0[1][1][0] };
// 	f32 p4[2] = { image1[1][0][0], image1[1][1][0] };

// 	f32 p5[2] = { image0[2][0][0], image0[2][1][0] };
// 	f32 p6[2] = { image1[2][0][0], image1[2][1][0] };

// 	Matrix K = Matrix(6, 9, {
//     p1[0], p1[1], 1, 0,0,0,-1*p1[0]*p2[0], -1*p1[1]*p2[0], -p2[0],
//     0,0,0,p1[0], p1[1], 1, -1*p1[0]*p2[1], -1*p1[1]*p2[1], -p2[1],
//     p3[0], p3[1], 1, 0,0,0,-1*p3[0]*p4[0], -1*p3[1]*p4[0], -p4[0],
//     0,0,0,p3[0], p3[1], 1, -1*p3[0]*p4[1], -1*p3[1]*p4[1], -p4[1],
//     p5[0], p5[1], 1, 0,0,0,-1*p5[0]*p6[0], -1*p5[1]*p6[0], -p6[0],
//     0,0,0,p5[0], p5[1], 1, -1*p5[0]*p6[1], -1*p5[1]*p6[1], -p6[1],
// 	});

// 	Matrix T = nullspace(K);

// 	Matrix U, S, V;

// 	printMatrix(T);
// 	svd(K, U, S, V);
// 	std::cout << "V:\n";
// 	printMatrix(V);
// 	std::cout << "S:\n";
// 	printMatrix(S);

// 	f32 B[9] = {
// 		nullSpaceVecB0(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecB1(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecB2(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecB3(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecB4(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecB5(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecB6(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecB7(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecB8(p1, p2, p3, p4, p5, p6)
// 	};

// 	f32 C[9] = {
// 		nullSpaceVecC0(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecC1(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecC2(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecC3(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecC4(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecC5(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecC6(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecC7(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecC8(p1, p2, p3, p4, p5, p6)
// 	};

// 	f32 D[9] = {
// 		nullSpaceVecD0(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecD1(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecD2(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecD3(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecD4(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecD5(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecD6(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecD7(p1, p2, p3, p4, p5, p6),
// 		nullSpaceVecD8(p1, p2, p3, p4, p5, p6)
// 	};

// 	for(i64 i=0; i<9; i++)
// 	{
// 		std::cout << B[i] << " ";
// 	}
// 	std::cout <<"\n";

// 	for(i64 i=0; i<9; i++)
// 	{
// 		std::cout << C[i] << " ";
// 	}
// 	std::cout <<"\n";

// 	for(i64 i=0; i<9; i++)
// 	{
// 		std::cout << D[i] << " ";
// 	}
// 	std::cout <<"\n";

// 	f32 dist = 10000000000;

// 	u64 best_i = 0;
// 	u64 best_j = 1;

// 	for(i64 i=0; i<3; i++)
// 	{
// 		for(i64 j=i; j<3; j++)
// 		{
// 			f32 k = norm(image0[i] - image0[j]);

// 			if(dist > k)
// 			{
// 				dist = k;
// 				best_i = i;
// 				best_j = j;
// 			}
// 		}
// 	}

// 	f32 pa1[2] = { image0[best_i][0][0], image0[best_i][1][0] };
// 	f32 pa2[2] = { image1[best_i][0][0], image1[best_i][1][0] };

// 	f32 pb1[2] = { image0[best_j][0][0], image0[best_j][1][0] };
// 	f32 pb2[2] = { image1[best_j][0][0], image1[best_j][1][0] };

// 	f32 a[2] = {angles[best_i], angles[best_j]};


// 	f32 b = betaP(B,C,D, pa1, pa2, pb1, pb2, a);
// 	f32 y = gammaP(B,C,D, pa1, pa2, pb1, pb2, a);

// 	std::cout << a[0] << "\n";
// 	std::cout << a[1] << "\n";
// 	std::cout << b << "\n";
// 	std::cout << y << "\n";

// 	f32 h0 = b*B[0] + y*C[0] + 1*D[0];
// 	f32 h1 = b*B[1] + y*C[1] + 1*D[1];
// 	f32 h2 = b*B[2] + y*C[2] + 1*D[2];
// 	f32 h3 = b*B[3] + y*C[3] + 1*D[3];
// 	f32 h4 = b*B[4] + y*C[4] + 1*D[4];
// 	f32 h5 = b*B[5] + y*C[5] + 1*D[5];
// 	f32 h6 = b*B[6] + y*C[6] + 1*D[6];
// 	f32 h7 = b*B[7] + y*C[7] + 1*D[7];
// 	f32 h8 = b*B[8] + y*C[8] + 1*D[8];

// 	return Matrix(3,3, {
// 		h0, h1, h2,
// 		h3, h4, h5,
// 		h6, h7, h8
// 	});
// }

f32 det3x3(Matrix& M)
{
	return 
		+ M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
		- M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
		+ M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
}

void fundamentalMatrixEstimation(Matrix& H, Matrix image0[2], Matrix image1[2])
{
	Matrix p1[7] = {
		Matrix(3,1, {10, 	 10, 1}),
		Matrix(3,1, {0, 	 10, 1}),
		Matrix(3,1, {10, 	 0,  1}),
		Matrix(3,1, {10, 	 5,  1}),
		Matrix(3,1, {5, 	 11, 1}),
		Matrix(3,1, {image0[1][0][0], 	 image0[1][1][0], 1}),
		Matrix(3,1, {image0[2][0][0], 	 image0[2][1][0], 1}),
	};

	Matrix p2[7] = {
		Matrix(0,0),
		Matrix(0,0),
		Matrix(0,0),
		Matrix(0,0),
		Matrix(0,0),
		Matrix(3,1, {image1[1][0][0], 	 image1[1][1][0], 1}),
		Matrix(3,1, {image1[2][0][0], 	 image1[2][1][0], 1}),
	};

	p2[0] = H*p1[0];
	p2[1] = H*p1[1];
	p2[2] = H*p1[2];
	p2[3] = H*p1[3];
	p2[4] = H*p1[4];

	Matrix D(7, 9, {
		p1[0][0][0]*p2[0][0][0], p1[0][1][0]*p2[0][0][0], p2[0][0][0], p1[0][0][0]*p2[0][1][0], p1[0][1][0]*p2[0][1][0], p2[0][1][0], p1[0][0][0], p1[0][1][0], 1,
		p1[1][0][0]*p2[1][0][0], p1[1][1][0]*p2[1][0][0], p2[1][0][0], p1[1][0][0]*p2[1][1][0], p1[1][1][0]*p2[1][1][0], p2[1][1][0], p1[1][0][0], p1[1][1][0], 1,
		p1[2][0][0]*p2[2][0][0], p1[2][1][0]*p2[2][0][0], p2[2][0][0], p1[2][0][0]*p2[2][1][0], p1[2][1][0]*p2[2][1][0], p2[2][1][0], p1[2][0][0], p1[2][1][0], 1,
		p1[3][0][0]*p2[3][0][0], p1[3][1][0]*p2[3][0][0], p2[3][0][0], p1[3][0][0]*p2[3][1][0], p1[3][1][0]*p2[3][1][0], p2[3][1][0], p1[3][0][0], p1[3][1][0], 1,
		p1[4][0][0]*p2[4][0][0], p1[4][1][0]*p2[4][0][0], p2[4][0][0], p1[4][0][0]*p2[4][1][0], p1[4][1][0]*p2[4][1][0], p2[4][1][0], p1[4][0][0], p1[4][1][0], 1,
		p1[5][0][0]*p2[5][0][0], p1[5][1][0]*p2[5][0][0], p2[5][0][0], p1[5][0][0]*p2[5][1][0], p1[5][1][0]*p2[5][1][0], p2[5][1][0], p1[5][0][0], p1[5][1][0], 1,
		p1[6][0][0]*p2[6][0][0], p1[6][1][0]*p2[6][0][0], p2[6][0][0], p1[6][0][0]*p2[6][1][0], p1[6][1][0]*p2[6][1][0], p2[6][1][0], p1[6][0][0], p1[6][1][0], 1,
	});

	// also possible with svd
	Matrix ns = nullspace(D);	
	Matrix f1(9,1, { ns[0][0], ns[0][1], ns[0][2], ns[0][3], ns[0][4], ns[0][5], ns[0][6], ns[0][7], ns[0][8] });
	Matrix f2(9,1, { ns[1][0], ns[1][1], ns[1][2], ns[1][3], ns[1][4], ns[1][5], ns[1][6], ns[1][7], ns[1][8] });

	Matrix F[2];

	F[0] = Matrix(3,3, {
		f1[0][0], f1[1][0], f1[2][0],
		f1[3][0], f1[4][0], f1[5][0],
		f1[6][0], f1[7][0], f1[8][0],
	});

	F[1] = Matrix(3,3, {
		f2[0][0], f2[1][0], f2[2][0],
		f2[3][0], f2[4][0], f2[5][0],
		f2[6][0], f2[7][0], f2[8][0],
	});

	f32 K[2][2][2];

	Matrix tmp(3,3);
	
	for(i64 i1=0; i1<2; i1++)
		for(i64 i2=0; i2<2; i2++)
			for(i64 i3=0; i3<2; i3++)
			{
				tmp[0][0] = F[i1][0][0];
				tmp[1][0] = F[i1][1][0];
				tmp[2][0] = F[i1][2][0];

				tmp[0][1] = F[i2][0][1];
				tmp[1][1] = F[i2][1][1];
				tmp[2][1] = F[i2][2][1];

				tmp[0][2] = F[i3][0][2];
				tmp[1][2] = F[i3][1][2];
				tmp[2][2] = F[i3][2][2];
		
				K[i1][i2][i3] = det3x3(tmp);
			}

	f32 coeffs[4];

	coeffs[0] = -K[1][0][0] + K[0][1][1] + K[0][0][0] + K[1][1][0] + K[1][0][1] - K[0][1][0] - K[0][0][1] - K[1][1][1];
	coeffs[1] = K[0][0][1] - 2*K[0][1][1] - 2*K[1][0][1]+K[1][0][0] - 2*K[1][1][0] + K[0][1][0] + 3*K[1][1][1];
	coeffs[2] = K[1][1][0] + K[0][1][1] + K[1][0][1] - 2*K[1][1][1];
	coeffs[3] = K[1][1][1];

	Polynomial poly(3, {coeffs[3], coeffs[2], coeffs[1], coeffs[0]});

	std::vector<f32> roots;
	
	poly.roots(roots);

	for(i64 i=0; i<roots.size(); i++)
	{
		Matrix Fundamental = F[0]*roots[i] + F[1]*(1-roots[i]);
		printMatrix(Fundamental);
	}
}


Matrix normalizePoints(Matrix* p, u64 n)
{
	f32 mx, my, sx, sy, s, msd;

	mx = 0.0;
	my = 0.0;

	for(int j=0; j<n; j++)
	{
		mx += p[j][0][0];
		my += p[j][1][0];
	}
	
	mx = mx/n;
	my = my/n;

	sx = 0.0;
	sy = 0.0;

	s = 0.0;

	for(int j=0; j<n; j++)
	{
		sx = pow(p[j][0][0] - mx, 2);
		sy = pow(p[j][1][0] - my, 2);
		s += sx + sy;
	}

	s = sqrt(s/(2 * n));

	Matrix T(3,3, {
		1./s, 0, -mx/s,
		0, 1./s, -my/s,
		0, 0, 1
	});

	for(int j=0; j<n; j++)
		p[j] = T*p[j];

	msd = 0.0;
	for(int j=0; j<n; j++)
		msd += pow(p[j][0][0],2) + pow(p[j][1][0],2);
	msd /= n;
	msd = sqrt(msd);

	// std::cout <<"msd^2: " << (msd*msd) << "\n";

	return T;
}

Matrix eightPointAlgorithm(Matrix x1[8], Matrix x2[8], bool normalized = true)
{
	Matrix lSingularVectors;
	Matrix singulaValues;
	Matrix rSingularVectors;

	Matrix A(8,9);
	Matrix T = identity(3,3);
	Matrix T_Inv = identity(3,3);

	if(normalized)
	{
		T 		= normalizePoints(x1, 8);
		T_Inv = normalizePoints(x2, 8);
	}

	// for(i32 i=0; i<8; i++)
	// {
	// 	A[i][0] = x1[i][0][0]*x2[i][0][0];
	// 	A[i][1] = x1[i][1][0]*x2[i][0][0];
	// 	A[i][2] = x2[i][0][0];
	// 	A[i][3] = x1[i][0][0]*x2[i][1][0];
	// 	A[i][4] = x1[i][1][0]*x2[i][1][0];
	// 	A[i][5] = x2[i][1][0];
	// 	A[i][6] = x1[i][0][0];
	// 	A[i][7] = x1[i][1][0];
	// 	A[i][8] = 1.0;
	// }

	for(i32 i=0; i<8; i++)
	{
		A[i][0] = x2[i][0][0]*x1[i][0][0];
		A[i][1] = x2[i][0][0]*x1[i][1][0];
		A[i][2] = x2[i][0][0]*x1[i][2][0];
		A[i][3] = x2[i][1][0]*x1[i][0][0];
		A[i][4] = x2[i][1][0]*x1[i][1][0];
		A[i][5] = x2[i][1][0]*x1[i][2][0];
		A[i][6] = x2[i][2][0]*x1[i][0][0];
		A[i][7] = x2[i][2][0]*x1[i][1][0];
		A[i][8] = x2[i][2][0]*x1[i][2][0];
	}

	Matrix U, V_T;

	f32* s = new f32[9];

	svd(A, U, s, V_T);

	Matrix F = Matrix(3,3);

	for(i32 i=0; i<3; i++)
	{
		for(i32 j=0; j<3; j++)
		{
			F[i][j] = V_T[8][i*3 + j]/V_T[8][8];
		}
	}

	svd(F, U, s, V_T);
	
	s[2] = 0;
	
	F = U*diag(s, 3, 3)*V_T;

	printMatrix(F);

	delete[] s;

	F = transpose(T_Inv)*F*T;
	F = F/F[2][2];

	return F;

}

Matrix getEssentialMatrix(Matrix F, Matrix K1, Matrix K2)
{
	return transpose(K2)*F*K1;
}

void estimateRotationAndTranslation(Matrix& E, Matrix& R1, Matrix& R2, Matrix& t1, Matrix& t2)
{
	f32 s[3]; 

	Matrix W = Matrix(3,3, {
		0, -1, 0,
		1, 0, 0,
		0, 0, 1
	});

	Matrix Z = Matrix(3,3, {
		0, 1, 0,
		-1, 0, 0,
		0, 0, 0
	});

	Matrix U, V_T;

	svd(E, U, s , V_T);

	f32 e = (s[0] + s[1]) / 2;
	
	s[0] = e;
	s[1] = e;
	s[2] = 0;

	E = U*diag(s, 3,3)*V_T;

	svd(E, U, s , V_T);

	R1 = U*W*V_T;
	R2 = U*transpose(W)*V_T;

	if(det3x3(R1) < 0)
	{
		R1 = R1 * -1;
	}

	if(det3x3(R2) < 0)
	{
		R2 = R2 * -1;
	}

	Matrix T = U*Z*transpose(U);

	Matrix t(3,1, {
		T[2][1],
		T[0][2],
		T[1][0],
	});

	t1 = t;

	t2 = t*-1;
}

f32 rad2deg(f32 rad)
{
	return rad * (180.0/3.141592653589793238463);
}

void rotationToEulerAngles(Matrix& R, f32* x, f32* y, f32* z)
{
	f32 sy = sqrt(R[0][0] * R[0][0] +  R[1][0] * R[1][0]);

  f32 singular = sy < 1e-6;

	if (!(sy < 1e-6))
	{
		*x = atan2(R[2][1] , R[2][2]);
		*y = atan2(-R[2][0], sy);
		*z = atan2(R[1][0], R[0][0]);
	}
	else
	{
		*x = atan2(-R[1][2], R[1][1]);
		*y = atan2(-R[2][0], sy);
		*z = 0;
	}

	std::cout << rad2deg(*x) << " " << rad2deg(*y) << " " << rad2deg(*z) << "\n";
}

Matrix triangulate(Matrix p, Matrix p_, Matrix P, Matrix P_)
{
	Matrix M(6,6);

	M[0][0] = P[0][0];
	M[1][0] = P[1][0];
	M[2][0] = P[2][0];

	M[0][1] = P[0][1];
	M[1][1] = P[1][1];
	M[2][1] = P[2][1];

	M[0][2] = P[0][2];
	M[1][2] = P[1][2];
	M[2][2] = P[2][2];

	M[0][3] = P[0][3];
	M[1][3] = P[1][3];
	M[2][3] = P[2][3];


	M[3][0] = P_[0][0];
	M[4][0] = P_[1][0];
	M[5][0] = P_[2][0];

	M[3][1] = P_[0][1];
	M[4][1] = P_[1][1];
	M[5][1] = P_[2][1];

	M[3][2] = P_[0][2];
	M[4][2] = P_[1][2];
	M[5][2] = P_[2][2];

	M[3][3] = P_[0][3];
	M[4][3] = P_[1][3];
	M[5][3] = P_[2][3];

	M[0][4] = -p[0][0];
	M[1][4] = -p[1][0];
	M[2][4] = -p[2][0];

	M[3][5] = -p_[0][0];
	M[4][5] = -p_[1][0];
	M[5][5] = -p_[2][0];

	Matrix U, V_T;

	f32 s[6];

	svd(M, U, s, V_T);

	Matrix X(4,1, {
		V_T[5][0]/V_T[5][3],
		V_T[5][1]/V_T[5][3],
		V_T[5][2]/V_T[5][3],
		V_T[5][3]/V_T[5][3]
	});

	// printMatrix(V_T);
	// std::cout << std::fixed << "(" << X[0][0] <<", " << X[1][0] << ", " << X[2][0] << ")\n";

	return X;
}

void chooseRealizableSolution(Matrix Intrinsics[2], Matrix Rotations[2], Matrix Translations[2], Matrix* pts1, Matrix* pts2, i32 n)
{
	i32 i, r, t, rotation, translation;

	Matrix cam0 = Intrinsics[0] * identity(3,4);
	Matrix cam1(3,4);
	
	bool realizable;

	for(r = 0; r < 2; r++)
	{
	
		realizable = true;
	
		for(t = 0; t < 2; t++)
		{
			std::cout << "TRIANGULATE: " << r << " " << t << "\n";
			
			cam1[0][0] = Rotations[r][0][0];
			cam1[1][0] = Rotations[r][1][0];
			cam1[2][0] = Rotations[r][2][0];
			
			cam1[0][1] = Rotations[r][0][1];
			cam1[1][1] = Rotations[r][1][1];
			cam1[2][1] = Rotations[r][2][1];

			cam1[0][2] = Rotations[r][0][2];
			cam1[1][2] = Rotations[r][1][2];
			cam1[2][2] = Rotations[r][2][2];

			cam1[0][3] = Translations[t][0][0];
			cam1[1][3] = Translations[t][1][0];
			cam1[2][3] = Translations[t][2][0];

			cam1 = Intrinsics[1] * cam1;

			for(i = 0; i < n; i++)
			{
				Matrix point = triangulate(pts1[i], pts2[i], cam0, cam1);

				// Reproject points		
				Matrix reprojected0 = cam0 * point;
				Matrix reprojected1 = cam1 * point;

				// if z of one of the reprojections are negative, the point
				// is behind the camera, so this configuration is invalid.
				// printMatrix(reprojected0);
				// std::cout << "\n";
				// printMatrix(reprojected1);
			
				if(reprojected0[2][0] < 0 || reprojected1[2][0] < 0)
				{
					realizable = false;
					break;
				}
				else
				{
					realizable = true;
				}
			}
		
			if(realizable)
			{
				rotation = r;
				translation = t;
			}
			if(!realizable)
			{
				break;
			}
		}
	

	}

	std::cout << rotation << " " << translation << "\n";
}


void estimateCameraFocalLengths(Matrix& Q)
{
	f32 s[4];

	Matrix U, W_T;
	svd(Q, U, s, W_T);

	if(det3x3(U) < 0)
	{
		U = U*-1;
	}

	if(det3x3(W_T) < 0)
	{
		W_T = W_T*-1;
	}

	Matrix E = Matrix(3,3, {
		0, 1, 0,
		-1, 0, 0,
		0, 0, 1
	});

	Matrix V = transpose(W_T)*E;

	std::cout << det3x3(U) << "\n";
	std::cout << det3x3(W_T) << "\n";

	double** U_ = new double*[3];
	U_[0] = new double[3];
	U_[1] = new double[3];
	U_[2] = new double[3];

	double** V_ = new double*[3];
	V_[0] = new double[3];
	V_[1] = new double[3];
	V_[2] = new double[3];

	U_[0][0] = U[0][0]; U_[0][1] = U[0][1]; U_[0][2] = U[0][2];
	U_[1][0] = U[1][0]; U_[1][1] = U[1][1]; U_[1][2] = U[1][2];
	U_[2][0] = U[2][0]; U_[2][1] = U[2][1]; U_[2][2] = U[2][2];

	V_[0][0] = V[0][0]; V_[0][1] = V[0][1]; V_[0][2] = V[0][2];
	V_[1][0] = V[1][0]; V_[1][1] = V[1][1]; V_[1][2] = V[1][2];
	V_[2][0] = V[2][0]; V_[2][1] = V[2][1]; V_[2][2] = V[2][2];

	double coeff0 = polyCoeff0(U_, V_, s[1], s[0]);
	double coeff1 = polyCoeff1(U_, V_, s[1], s[0]);
	double coeff2 = polyCoeff2(U_, V_, s[1], s[0]);
	double coeff3 = polyCoeff3(U_, V_, s[1], s[0]);

	printf("%f, %f, %f, %f\n", coeff0, coeff1, coeff2, coeff3);
	std::cout << std::scientific << coeff0 << " " << coeff1 << " " << coeff2 << " " << coeff3 << std::endl; 

	delete[] U_[0];
	delete[] U_[1];
	delete[] U_[2];
	delete[] U_;

	delete[] V_[0];
	delete[] V_[1];
	delete[] V_[2];
	delete[] V_;
}
