#pragma once

#include <cmath>
#include "bundle/codegen/Homography.hpp"
#include "algebra/matrix/Matrix.hpp"
#include "algebra/polynomial/Polynomial.hpp"
#include "algebra/linear/Linear.hpp"
#include "algebra/linear/SingularValueDecomposition.hpp"
#include "algebra/linear/svd.hpp"

#include "algebra/linear/Eigen/SVD"

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

Matrix coplanarThreePointHomography(Matrix image0[3], Matrix image1[3], f32 angles[3])
{
	f32 p1[2] = { image0[0][0][0], image0[0][1][0] };
	f32 p2[2] = { image1[0][0][0], image1[0][1][0] };

	f32 p3[2] = { image0[1][0][0], image0[1][1][0] };
	f32 p4[2] = { image1[1][0][0], image1[1][1][0] };

	f32 p5[2] = { image0[2][0][0], image0[2][1][0] };
	f32 p6[2] = { image1[2][0][0], image1[2][1][0] };

	Matrix K = Matrix(6, 9, {
    p1[0], p1[1], 1, 0,0,0,-1*p1[0]*p2[0], -1*p1[1]*p2[0], -p2[0],
    0,0,0,p1[0], p1[1], 1, -1*p1[0]*p2[1], -1*p1[1]*p2[1], -p2[1],
    p3[0], p3[1], 1, 0,0,0,-1*p3[0]*p4[0], -1*p3[1]*p4[0], -p4[0],
    0,0,0,p3[0], p3[1], 1, -1*p3[0]*p4[1], -1*p3[1]*p4[1], -p4[1],
    p5[0], p5[1], 1, 0,0,0,-1*p5[0]*p6[0], -1*p5[1]*p6[0], -p6[0],
    0,0,0,p5[0], p5[1], 1, -1*p5[0]*p6[1], -1*p5[1]*p6[1], -p6[1],
	});

	Matrix T = nullspace(K);

	Matrix U, S, V;

	printMatrix(T);
	svd(K, U, S, V);
	std::cout << "V:\n";
	printMatrix(V);
	std::cout << "S:\n";
	printMatrix(S);

	f32 B[9] = {
		nullSpaceVecB0(p1, p2, p3, p4, p5, p6),
		nullSpaceVecB1(p1, p2, p3, p4, p5, p6),
		nullSpaceVecB2(p1, p2, p3, p4, p5, p6),
		nullSpaceVecB3(p1, p2, p3, p4, p5, p6),
		nullSpaceVecB4(p1, p2, p3, p4, p5, p6),
		nullSpaceVecB5(p1, p2, p3, p4, p5, p6),
		nullSpaceVecB6(p1, p2, p3, p4, p5, p6),
		nullSpaceVecB7(p1, p2, p3, p4, p5, p6),
		nullSpaceVecB8(p1, p2, p3, p4, p5, p6)
	};

	f32 C[9] = {
		nullSpaceVecC0(p1, p2, p3, p4, p5, p6),
		nullSpaceVecC1(p1, p2, p3, p4, p5, p6),
		nullSpaceVecC2(p1, p2, p3, p4, p5, p6),
		nullSpaceVecC3(p1, p2, p3, p4, p5, p6),
		nullSpaceVecC4(p1, p2, p3, p4, p5, p6),
		nullSpaceVecC5(p1, p2, p3, p4, p5, p6),
		nullSpaceVecC6(p1, p2, p3, p4, p5, p6),
		nullSpaceVecC7(p1, p2, p3, p4, p5, p6),
		nullSpaceVecC8(p1, p2, p3, p4, p5, p6)
	};

	f32 D[9] = {
		nullSpaceVecD0(p1, p2, p3, p4, p5, p6),
		nullSpaceVecD1(p1, p2, p3, p4, p5, p6),
		nullSpaceVecD2(p1, p2, p3, p4, p5, p6),
		nullSpaceVecD3(p1, p2, p3, p4, p5, p6),
		nullSpaceVecD4(p1, p2, p3, p4, p5, p6),
		nullSpaceVecD5(p1, p2, p3, p4, p5, p6),
		nullSpaceVecD6(p1, p2, p3, p4, p5, p6),
		nullSpaceVecD7(p1, p2, p3, p4, p5, p6),
		nullSpaceVecD8(p1, p2, p3, p4, p5, p6)
	};

	for(i64 i=0; i<9; i++)
	{
		std::cout << B[i] << " ";
	}
	std::cout <<"\n";

	for(i64 i=0; i<9; i++)
	{
		std::cout << C[i] << " ";
	}
	std::cout <<"\n";

	for(i64 i=0; i<9; i++)
	{
		std::cout << D[i] << " ";
	}
	std::cout <<"\n";

	f32 dist = 10000000000;

	u64 best_i = 0;
	u64 best_j = 1;

	for(i64 i=0; i<3; i++)
	{
		for(i64 j=i; j<3; j++)
		{
			f32 k = norm(image0[i] - image0[j]);

			if(dist > k)
			{
				dist = k;
				best_i = i;
				best_j = j;
			}
		}
	}

	f32 pa1[2] = { image0[best_i][0][0], image0[best_i][1][0] };
	f32 pa2[2] = { image1[best_i][0][0], image1[best_i][1][0] };

	f32 pb1[2] = { image0[best_j][0][0], image0[best_j][1][0] };
	f32 pb2[2] = { image1[best_j][0][0], image1[best_j][1][0] };

	f32 a[2] = {angles[best_i], angles[best_j]};


	f32 b = betaP(B,C,D, pa1, pa2, pb1, pb2, a);
	f32 y = gammaP(B,C,D, pa1, pa2, pb1, pb2, a);

	std::cout << a[0] << "\n";
	std::cout << a[1] << "\n";
	std::cout << b << "\n";
	std::cout << y << "\n";

	f32 h0 = b*B[0] + y*C[0] + 1*D[0];
	f32 h1 = b*B[1] + y*C[1] + 1*D[1];
	f32 h2 = b*B[2] + y*C[2] + 1*D[2];
	f32 h3 = b*B[3] + y*C[3] + 1*D[3];
	f32 h4 = b*B[4] + y*C[4] + 1*D[4];
	f32 h5 = b*B[5] + y*C[5] + 1*D[5];
	f32 h6 = b*B[6] + y*C[6] + 1*D[6];
	f32 h7 = b*B[7] + y*C[7] + 1*D[7];
	f32 h8 = b*B[8] + y*C[8] + 1*D[8];

	return Matrix(3,3, {
		h0, h1, h2,
		h3, h4, h5,
		h6, h7, h8
	});
}

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
	f32 tx, ty, sx, sy, s, msd;

	tx = 0.0;
	ty = 0.0;

	for(int j=0; j<n; j++)
	{
		tx += p[j][0][0];
		ty += p[j][1][0];
	}
	
	tx = tx/(n);
	ty = ty/(n);

	// for(int j=0; j<n; j++)
	// {
	// 	p[j][0][0] = (p[j][0][0] - tx);
	// 	p[j][1][0] = (p[j][1][0] - ty);
	// }

	sx = 0.0;
	sy = 0.0;

	s = 0.0;

	for(int j=0; j<n; j++)
	{
		sx = pow(p[j][0][0] - tx,2);
		sy = pow(p[j][1][0] - ty,2);
		s += sx+sy;
	}

	s /= 2*n;
	s = sqrt(s);

	Matrix T(3,3, {
		1./s, 0, -tx/s,
		0, 1./s, -ty/s,
		0, 0, 1
	});

	for(int j=0; j<n; j++)
		p[j] = T*p[j];

	msd = 0.0;
	for(int j=0; j<n; j++)
		msd += pow(p[j][0][0],2) + pow(p[j][1][0],2);
	msd /= n;
	msd = sqrt(msd);

	std::cout <<"msd^2: " << (msd*msd) << "\n";

	return T;
}

Matrix eightPointAlgorithm(Matrix image0[8], Matrix image1[8])
{

	Matrix lSingularVectors;
	Matrix singulaValues;
	Matrix rSingularVectors;

	Matrix p[8] = {
		Matrix(3,1, {image0[0][0][0], 	 image0[0][1][0], 1}),
		Matrix(3,1, {image0[1][0][0], 	 image0[1][1][0], 1}),
		Matrix(3,1, {image0[2][0][0], 	 image0[2][1][0], 1}),
		Matrix(3,1, {image0[3][0][0], 	 image0[3][1][0], 1}),
		Matrix(3,1, {image0[4][0][0], 	 image0[4][1][0], 1}),
		Matrix(3,1, {image0[5][0][0], 	 image0[5][1][0], 1}),
		Matrix(3,1, {image0[6][0][0], 	 image0[6][1][0], 1}),
		Matrix(3,1, {image0[7][0][0], 	 image0[7][1][0], 1}),
	};

	Matrix p_[8] = {
		Matrix(3,1, {image1[0][0][0], 	 image1[0][1][0], 1}),
		Matrix(3,1, {image1[1][0][0], 	 image1[1][1][0], 1}),
		Matrix(3,1, {image1[2][0][0], 	 image1[2][1][0], 1}),
		Matrix(3,1, {image1[3][0][0], 	 image1[3][1][0], 1}),
		Matrix(3,1, {image1[4][0][0], 	 image1[4][1][0], 1}),
		Matrix(3,1, {image1[5][0][0], 	 image1[5][1][0], 1}),
		Matrix(3,1, {image1[6][0][0], 	 image1[6][1][0], 1}),
		Matrix(3,1, {image1[7][0][0], 	 image1[7][1][0], 1}),
	};

	// normalize points
	// Matrix Tinv = normalizePoints(p, 8);
	// Matrix T_inv = normalizePoints(p_, 8);

	// Matrix A(8, 9, {
	// 	p[0][0][0]*p_[0][0][0], p[0][1][0]*p_[0][0][0], p_[0][0][0], p[0][0][0]*p_[0][1][0], p[0][1][0]*p_[0][1][0], p_[0][1][0], p[0][0][0],  p[0][1][0], 1,
	// 	p[1][0][0]*p_[1][0][0], p[1][1][0]*p_[1][0][0], p_[1][0][0], p[1][0][0]*p_[1][1][0], p[1][1][0]*p_[1][1][0], p_[1][1][0], p[1][0][0],  p[1][1][0], 1,
	// 	p[2][0][0]*p_[2][0][0], p[2][1][0]*p_[2][0][0], p_[2][0][0], p[2][0][0]*p_[2][1][0], p[2][1][0]*p_[2][1][0], p_[2][1][0], p[2][0][0],  p[2][1][0], 1,
	// 	p[3][0][0]*p_[3][0][0], p[3][1][0]*p_[3][0][0], p_[3][0][0], p[3][0][0]*p_[3][1][0], p[3][1][0]*p_[3][1][0], p_[3][1][0], p[3][0][0],  p[3][1][0], 1,
	// 	p[4][0][0]*p_[4][0][0], p[4][1][0]*p_[4][0][0], p_[4][0][0], p[4][0][0]*p_[4][1][0], p[4][1][0]*p_[4][1][0], p_[4][1][0], p[4][0][0],  p[4][1][0], 1,
	// 	p[5][0][0]*p_[5][0][0], p[5][1][0]*p_[5][0][0], p_[5][0][0], p[5][0][0]*p_[5][1][0], p[5][1][0]*p_[5][1][0], p_[5][1][0], p[5][0][0],  p[5][1][0], 1,
	// 	p[6][0][0]*p_[6][0][0], p[6][1][0]*p_[6][0][0], p_[6][0][0], p[6][0][0]*p_[6][1][0], p[6][1][0]*p_[6][1][0], p_[6][1][0], p[6][0][0],  p[6][1][0], 1,
	// 	p[7][0][0]*p_[7][0][0], p[7][1][0]*p_[7][0][0], p_[7][0][0], p[7][0][0]*p_[7][1][0], p[7][1][0]*p_[7][1][0], p_[7][1][0], p[7][0][0],  p[7][1][0], 1,
	// });
	
	Eigen::MatrixXf A(8,9);

	A << p_[0][0][0]*p[0][0][0], p_[0][0][0]*p[0][1][0], p_[0][0][0], p[0][0][0]*p_[0][1][0], p[0][1][0]*p_[0][1][0], p_[0][1][0], p[0][0][0],  p[0][1][0], 1,
	p_[1][0][0]*p[1][0][0], p_[1][0][0]*p[1][1][0], p_[1][0][0], p[1][0][0]*p_[1][1][0], p[1][1][0]*p_[1][1][0], p_[1][1][0], p[1][0][0],  p[1][1][0], 1,
	p_[2][0][0]*p[2][0][0], p_[2][0][0]*p[2][1][0], p_[2][0][0], p[2][0][0]*p_[2][1][0], p[2][1][0]*p_[2][1][0], p_[2][1][0], p[2][0][0],  p[2][1][0], 1,
	p_[3][0][0]*p[3][0][0], p_[3][0][0]*p[3][1][0], p_[3][0][0], p[3][0][0]*p_[3][1][0], p[3][1][0]*p_[3][1][0], p_[3][1][0], p[3][0][0],  p[3][1][0], 1,
	p_[4][0][0]*p[4][0][0], p_[4][0][0]*p[4][1][0], p_[4][0][0], p[4][0][0]*p_[4][1][0], p[4][1][0]*p_[4][1][0], p_[4][1][0], p[4][0][0],  p[4][1][0], 1,
	p_[5][0][0]*p[5][0][0], p_[5][0][0]*p[5][1][0], p_[5][0][0], p[5][0][0]*p_[5][1][0], p[5][1][0]*p_[5][1][0], p_[5][1][0], p[5][0][0],  p[5][1][0], 1,
	p_[6][0][0]*p[6][0][0], p_[6][0][0]*p[6][1][0], p_[6][0][0], p[6][0][0]*p_[6][1][0], p[6][1][0]*p_[6][1][0], p_[6][1][0], p[6][0][0],  p[6][1][0], 1,
	p_[7][0][0]*p[7][0][0], p_[7][0][0]*p[7][1][0], p_[7][0][0], p[7][0][0]*p_[7][1][0], p[7][1][0]*p_[7][1][0], p_[7][1][0], p[7][0][0],  p[7][1][0], 1;

	// A << p[0][0][0]*p_[0][0][0], p[0][1][0]*p_[0][0][0], p_[0][0][0], p[0][0][0]*p_[0][1][0], p[0][1][0]*p_[0][1][0], p_[0][1][0], p[0][0][0],  p[0][1][0], 1,
	// p[1][0][0]*p_[1][0][0], p[1][1][0]*p_[1][0][0], p_[1][0][0], p[1][0][0]*p_[1][1][0], p[1][1][0]*p_[1][1][0], p_[1][1][0], p[1][0][0],  p[1][1][0], 1,
	// p[2][0][0]*p_[2][0][0], p[2][1][0]*p_[2][0][0], p_[2][0][0], p[2][0][0]*p_[2][1][0], p[2][1][0]*p_[2][1][0], p_[2][1][0], p[2][0][0],  p[2][1][0], 1,
	// p[3][0][0]*p_[3][0][0], p[3][1][0]*p_[3][0][0], p_[3][0][0], p[3][0][0]*p_[3][1][0], p[3][1][0]*p_[3][1][0], p_[3][1][0], p[3][0][0],  p[3][1][0], 1,
	// p[4][0][0]*p_[4][0][0], p[4][1][0]*p_[4][0][0], p_[4][0][0], p[4][0][0]*p_[4][1][0], p[4][1][0]*p_[4][1][0], p_[4][1][0], p[4][0][0],  p[4][1][0], 1,
	// p[5][0][0]*p_[5][0][0], p[5][1][0]*p_[5][0][0], p_[5][0][0], p[5][0][0]*p_[5][1][0], p[5][1][0]*p_[5][1][0], p_[5][1][0], p[5][0][0],  p[5][1][0], 1,
	// p[6][0][0]*p_[6][0][0], p[6][1][0]*p_[6][0][0], p_[6][0][0], p[6][0][0]*p_[6][1][0], p[6][1][0]*p_[6][1][0], p_[6][1][0], p[6][0][0],  p[6][1][0], 1,
	// p[7][0][0]*p_[7][0][0], p[7][1][0]*p_[7][0][0], p_[7][0][0], p[7][0][0]*p_[7][1][0], p[7][1][0]*p_[7][1][0], p_[7][1][0], p[7][0][0],  p[7][1][0], 1;

	Eigen::JacobiSVD<Eigen::MatrixXf> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

	Matrix Fnorm = Matrix(3,3, {
		svd.matrixV()(0,8), svd.matrixV()(1,8), svd.matrixV()(2,8),
		svd.matrixV()(3,8), svd.matrixV()(4,8), svd.matrixV()(5,8),
		svd.matrixV()(6,8), svd.matrixV()(7,8), svd.matrixV()(8,8),
	});

	// Eigen::MatrixXf U = svd.matrixU();
	// Eigen::MatrixXf D = svd.singularValues().matrix().asDiagonal();
	// Eigen::MatrixXf V = svd.matrixV();

	// std::cout << A << "\n";
	// std::cout << U << "\n";
	// std::cout << D << "\n";
	// std::cout << V << "\n";
	// std::cout << a_ - (U*D) << "\n";
	// Matrix Fnorm = Matrix(3,3, {
	// 	rSingularVectors[0][8], rSingularVectors[1][8], rSingularVectors[2][8],
	// 	rSingularVectors[3][8], rSingularVectors[4][8], rSingularVectors[5][8],
	// 	rSingularVectors[6][8], rSingularVectors[7][8], rSingularVectors[8][8],
	// });

	// jacobiSVD(Fnorm, singulaValues, lSingularVectors, rSingularVectors, 100);

	// Matrix S(3,3, {
	// 	singulaValues[0][0], 0, 0,
	// 	0, singulaValues[1][0], 0,
	// 	0, 0, singulaValues[2][0]
	// });

 	// Matrix _Fnorm = lSingularVectors*S*transpose(rSingularVectors);

	// printMatrix(Fnorm);
	// std::cout << "ASDASDSADAS\n";
	// printMatrix(Fnorm - _Fnorm);
	// std::cout << "ASDASDSADAS\n";

	// printMatrix(Fnorm);

	// for(i64 i=0; i<8; i++)
	// {
	// 	Matrix err = transpose(p[i])*Fnorm*p_[i];
	// 	printMatrix(err);
	// }


	return Fnorm;
	// Matrix F = transpose(Tinv)*Fnorm*T_inv;

	// std::cout <<"det: " <<  det3x3(Fnorm) << "\n";
	// std::cout <<"det: " <<  det3x3(F) << "\n";

	// return F;
	// return F;
	// Matrix S,V,D;
	// svd(A, S, D, V);
	// printMatrix(A);
	// std::cout <<"ASDASDS\n";
	// printMatrix(S);
	// std::cout <<"ASDASDS\n";
	// printMatrix(D);
	// std::cout <<"ASDASDS\n";
	// printMatrix(V);
	// std::cout <<"ASDASDS\n";
	// printMatrix(ns);
	// std::cout <<"*************\n";

	// Matrix W(9,9, {
	// 	D[0][0], 0, 0, 0, 0, 0, 0, 0, 0,
	// 	0, D[1][0], 0, 0, 0, 0, 0, 0, 0,
	// 	0, 0, D[2][0], 0, 0, 0, 0, 0, 0,
	// 	0, 0, 0, D[3][0], 0, 0, 0, 0, 0,
	// 	0, 0, 0, 0, D[4][0], 0, 0, 0, 0,
	// 	0, 0, 0, 0, 0, D[5][0], 0, 0, 0,
	// 	0, 0, 0, 0, 0, 0, D[6][0], 0, 0,
	// 	0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 	0, 0, 0, 0, 0, 0, 0, 0, 0,
	// });
	// Matrix V_T = transpose(V);
	// Matrix A_ = S*W*V_T;

	// printMatrix(A_);

	// Matrix f1(9,1, { ns[0][0], ns[0][1], ns[0][2], ns[0][3], ns[0][4], ns[0][5], ns[0][6], ns[0][7], ns[0][8] });
	// Matrix f2(9,1, { ns[1][0], ns[1][1], ns[1][2], ns[1][3], ns[1][4], ns[1][5], ns[1][6], ns[1][7], ns[1][8] });

	// Matrix F[2];

	// F[0] = Matrix(3,3, {
	// 	f1[0][0], f1[1][0], f1[2][0],
	// 	f1[3][0], f1[4][0], f1[5][0],
	// 	f1[6][0], f1[7][0], f1[8][0],
	// });

	// F[1] = Matrix(3,3, {
	// 	f2[0][0], f2[1][0], f2[2][0],
	// 	f2[3][0], f2[4][0], f2[5][0],
	// 	f2[6][0], f2[7][0], f2[8][0],
	// });

	// f32 K[2][2][2];

	// Matrix tmp(3,3);
	
	// for(i64 i1=0; i1<2; i1++)
	// 	for(i64 i2=0; i2<2; i2++)
	// 		for(i64 i3=0; i3<2; i3++)
	// 		{
	// 			tmp[0][0] = F[i1][0][0];
	// 			tmp[1][0] = F[i1][1][0];
	// 			tmp[2][0] = F[i1][2][0];

	// 			tmp[0][1] = F[i2][0][1];
	// 			tmp[1][1] = F[i2][1][1];
	// 			tmp[2][1] = F[i2][2][1];

	// 			tmp[0][2] = F[i3][0][2];
	// 			tmp[1][2] = F[i3][1][2];
	// 			tmp[2][2] = F[i3][2][2];
		
	// 			K[i1][i2][i3] = det3x3(tmp);
	// 		}

	// f32 coeffs[4];

	// coeffs[0] = -K[1][0][0] + K[0][1][1] + K[0][0][0] + K[1][1][0] + K[1][0][1] - K[0][1][0] - K[0][0][1] - K[1][1][1];
	// coeffs[1] = K[0][0][1] - 2*K[0][1][1] - 2*K[1][0][1]+K[1][0][0] - 2*K[1][1][0] + K[0][1][0] + 3*K[1][1][1];
	// coeffs[2] = K[1][1][0] + K[0][1][1] + K[1][0][1] - 2*K[1][1][1];
	// coeffs[3] = K[1][1][1];

	// Polynomial poly(3, {coeffs[3], coeffs[2], coeffs[1], coeffs[0]});

	// std::vector<f32> roots;
	
	// poly.roots(roots, 0.00000000000001);

	// // printPoly(poly);

	// for(i64 i=0; i<roots.size(); i++)
	// {
	// 	std::cout << "ASDASD\n";
	// 	Matrix Fundamental = F[0]*roots[i] + F[1]*(1-roots[i]);
	
	// 	printMatrix(Fundamental);
	
	// 	// Matrix U,W,V;
	// 	// svd(Fundamental, U, W, V);

	// 	// Matrix D(3,3, {
	// 	// 	W[0][0], 0, 0,
	// 	// 	0, W[1][0], 0,
	// 	// 	0, 0, W[2][0],
	// 	// });
	// 	// if(det3x3(U) == -1.)
	// 	// 	U = U*-1;
	// 	// if(det3x3(V) == -1.)
	// 	// 	V = V*-1;
	// 	// std::cout << det3x3(U) << "\n";
	// 	// std::cout << det3x3(V) << "\n";
	// 	// printMatrix(U);
	// 	// std::cout << "\n";
	// 	// printMatrix(W);
	// 	// std::cout << "\n";
	// 	// printMatrix(V);

	// 	// Matrix E(3,3, {
	// 	// 	0, 1, 0,
	// 	// 	-1, 0, 0,
	// 	// 	0, 0, 1
	// 	// });

	// 	// Matrix Z(3,3, {
	// 	// 	0, -1, 0,
	// 	// 	1, 0, 0,
	// 	// 	0, 0, 0
	// 	// });

	// 	// Matrix X0 = U*E*transpose(V);
	// 	// Matrix X1 = U*transpose(E)*transpose(V);
	
	// 	// Matrix v(3,1, {0,0,1});
	
	// 	// Matrix t0 = U*v;
	// 	// Matrix t1 = U*-1*v;
	// 	// t0 = transpose(t0);
	// 	// t1 = transpose(t1);
	// 	// printMatrix(t0);
	// 	// std::cout <<"\n";
	// 	// printMatrix(X0);
	// 	// std::cout <<"\n";
	// 	// std::cout <<"\n";
	// 	// printMatrix(t1);
	// 	// std::cout <<"\n";
	// 	// printMatrix(X0);
	// 	// std::cout <<"\n";
	// 	// std::cout <<"\n";
	// 	// printMatrix(t0);
	// 	// std::cout <<"\n";
	// 	// printMatrix(X1);
	// 	// std::cout <<"\n";
	// 	// std::cout <<"\n";
	// 	// printMatrix(t1);
	// 	// std::cout <<"\n";
	// 	// printMatrix(X1);
	// 	for(i64 j=0; j<7; j++)
	// 	{
	// 		Matrix k = transpose(p_[j])*Fundamental*p[j];
	// 		printMatrix(k);
	// 	}
	
	// 	std::cout << "\n";
	// }
}

