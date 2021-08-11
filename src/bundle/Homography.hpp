#pragma once

#include <assert.h>

#include "algebra/core/types.hpp"
#include "algebra/matrix/Matrix.hpp"
#include "algebra/linear/Linear.hpp"
#include "bundle/codegen/Homography.hpp"

using namespace karu;
using namespace karu::algebra;

Matrix coplanarThreePointHomography(Matrix image0[3], Matrix image1[3], f32 angles[3])
{
	f32 p1[2] = { image0[0][0][0], image0[0][1][0] };
	f32 p2[2] = { image1[0][0][0], image1[0][1][0] };

	f32 p3[2] = { image0[1][0][0], image0[1][1][0] };
	f32 p4[2] = { image1[1][0][0], image1[1][1][0] };

	f32 p5[2] = { image0[2][0][0], image0[2][1][0] };
	f32 p6[2] = { image1[2][0][0], image1[2][1][0] };

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

	f32 dist = 340282350000000000000000000000000000000;

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

	f32 b = beta(B,C,D, pa1, pa2, pb1, pb2, a);
	f32 y = gamma(B,C,D, pa1, pa2, pb1, pb2, a);

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

	Matrix Dns = nullspace(D);
	
	Matrix e(9,1, {Dns[0][0], Dns[0][1], Dns[0][2], Dns[0][3], Dns[0][4], Dns[0][5], Dns[0][6], Dns[0][7], Dns[0][8]});
	Matrix g(9,1, {Dns[1][0], Dns[1][1], Dns[1][2], Dns[1][3], Dns[1][4], Dns[1][5], Dns[1][6], Dns[1][7], Dns[1][8]});

	F = x*e + (1-x)*g;
}
