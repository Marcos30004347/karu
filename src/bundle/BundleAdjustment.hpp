/**
** Bundle Adjustment Basics:

	problem:
		min f 
		where f = sum of i,m {(u[i][j] - proj(C[j], X[i]))**2}

		where u[i][j] is the observed point X[i] by the
		camera C[j].

	def:
		r = u[i][j] - proj(C[j], X[i])
	
		given x = [c, p] with 
		c being a flat vector with all the camera parameters and 
		p being a flat vector with all the point parameters. 

		lets define û[i][j] = proj(C[j], X[i])
		then define:
		 	A[i][j] = diff(û[i][j], c[j]), and diff(û[i][j], c[k]) = 0 for every j different from k 
 		 	B[i][j] = diff(û[i][j], p[i]), and diff(û[i][j], p[k]) = 0 for every i different from k 

	U[j] = sum over i for each point i that is seen by camera j of 	A[i][j].T * A[i][j]
	V[i] = sum over j for each camera j that sees the point i of 		B[i][j].T * B[i][j]
	W[i][j] = A[i][j].T * B[i][j]

	r[i][j] = u[i][j] - û[i][j]
	e[c[j]] = sum over all i seen in camera j of A[i][j]*r[i][j]
	e[p[i]] = sum over all j that sees point i of B[i][j]*r[i][j]

	Then we can get delta x by
	[[U W], [W.T, V]] * delta(x) = flat([e[c], e[p]])

	It is possibly to get 
	(U - W*inverse(V)*W.T) = S
	with being c is 1: Compute the derivation matrix:
Aij =
∂uˆij
∂cj
, Bij =
∂uˆij
∂Xi
rij = uij − π(cj , Xi) = uij − uˆij


/**
Algorithm: Compute Sx with W, V and x
without forming S on memory:
	x1 = W.T*x
	x2 = inverse(V) * x1
	x3 = W*x2
	x4 = B*x
	Sx = x4 - x3
**/
#pragma once

#include <cmath>
#include <vector>
#include <set>
#include <array>
#include <assert.h>
#include "algebra/matrix/Matrix.hpp"
#include "algebra/sparse/SpMatrix.hpp"
#include "camera/Camera.hpp"

using namespace karu::algebra;

namespace karu::bundle 
{

struct Point 
{
	f32 x;
	f32 y;
	f32 z;
};

struct Pixel 
{
	f32 u;
	f32 v;
};

struct Bundle 
{
	// Camera parameter
	Camera camera;

	// Camera observerd points pixel coordinates
	std::vector<Pixel> projections;

	// idx of the point that projections[i] is the projection in this camera
	std::vector<u64> point_idx;
};

Matrix packObservations(std::vector<Bundle>& bundles, std::vector<Point>& points)
{
	std::vector<f32> x(0);
	x.reserve(15*bundles.size() + 3*points.size());

	for(Bundle b : bundles)
	{
		x.push_back(b.camera.fx);
		x.push_back(b.camera.fy);
		x.push_back(b.camera.cx);
		x.push_back(b.camera.cy);
		x.push_back(b.camera.P[0][0]);
		x.push_back(b.camera.P[1][0]);
		x.push_back(b.camera.P[2][0]);
		x.push_back(b.camera.R[0][0]);
		x.push_back(b.camera.R[1][0]);
		x.push_back(b.camera.R[2][0]);
		x.push_back(b.camera.k1);
		x.push_back(b.camera.k2);
		x.push_back(b.camera.k3);
		x.push_back(b.camera.p1);
		x.push_back(b.camera.p2);
	}

	for(Point p : points)
	{
		x.push_back(p.x);
		x.push_back(p.y);
		x.push_back(p.z);
	}

	return Matrix(15, 1, x.data());
}

// Point i, camera j
Matrix A(std::vector<Bundle>& bundles, std::vector<Matrix>& points, u64 i, u64 j)
{
	return Matrix(2, 15, {
		u_dfx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dfy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dcx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dcy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dCx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dCy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dCz(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dr1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dr2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dr3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dk1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dk2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dk3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dp1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dp2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dfx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dfy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dcx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dcy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dCx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dCy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dCz(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dr1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dr2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dr3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dk1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dk2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dk3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dp1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dp2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
	});
}

Matrix B(std::vector<Bundle>& bundles, std::vector<Matrix>& points, u64 i, u64 j)
{
	return Matrix(2, 3, {
		u_dX(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dY(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		u_dZ(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),

		v_dX(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dY(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		v_dZ(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
	});
}

void build_U(std::vector<Bundle>& bundles, std::vector<Matrix>& points, SpMatrix& U)
{
	std::vector<u64> U_rows;
	std::vector<u64> U_cols_idx;
	
	for(i64 j=0; j<bundles.size(); j++)
	{
		U_rows.push_back(j);
		U_cols_idx.push_back(j*15);
	}

	U_rows.push_back(bundles.size());

	std::vector<f32> U_data(15*15*U_cols_idx.size(), 0);

	for(i64 j=0; j<bundles.size(); j++)
	{
		for(u64 i : bundles[j].point_idx)
		{
			Matrix A_block = A(bundles, points, i, j);
	
			Matrix U_block = transpose(A_block) * A_block; 

			for(i64 c=0; c<15; c++)
			{
				for(i64 l=0; l<15; l++)
				{
					U_data[j*15*15 + l + c*15] += U_block[l][c];
				}
			}
		}
	}

	U = SpMatrix((U_rows.size() - 1)*15, 15*U_cols_idx.size(), 15, 15, U_rows, U_cols_idx, U_data);
}

void build_V(std::vector<Bundle>& bundles, std::vector<Matrix>& points, SpMatrix& V)
{
	std::vector<u64> V_rows;
	std::vector<u64> V_cols_idx;

	for(i64 j=0; j<points.size(); j++)
	{
		V_rows.push_back(j);
		V_cols_idx.push_back(j*3);
	}

	V_rows.push_back(points.size());

	std::vector<f32> V_data(9*V_cols_idx.size(), 0);
	
	for(i64 j=0; j<bundles.size(); j++)
	{
		for(u64 i : bundles[j].point_idx)
		{
			Matrix B_block = B(bundles, points, i, j);
			Matrix V_block = transpose(B_block) * B_block; 

			// Push data into col major layout
			for(i64 c = 0; c<3; c++)
			{
				for(i64 l = 0; l<3; l++)
				{
					V_data[i*9 + l + c*3] += V_block[l][c];
				}
			}
		}
	}
	V = SpMatrix((V_rows.size() - 1)*3, 3*V_cols_idx.size(), 3, 3, V_rows, V_cols_idx, V_data);
}

void build_W(std::vector<Bundle>& bundles, std::vector<Matrix>& points, SpMatrix& W)
{
	u64 W_block_cnt = 0;

	for(i64 i=0; i<bundles.size(); i++)
		for(u64 j : bundles[i].point_idx)
			W_block_cnt++;

	std::vector<f32> W_data(45*W_block_cnt, 0);
	std::vector<u64> W_cols_idx(W_block_cnt, 0);
	std::vector<u64> W_rows(bundles.size()+1, 0);

	u64 idx = 0;

	for(i64 i=0; i<bundles.size(); i++)
	{
		for(u64 j : bundles[i].point_idx)
		{
			W_cols_idx[idx++] = 3*j;
			W_rows[i+1]++;
		}

		if(i+2 < W_rows.size())
			W_rows[i+2] = W_rows[i+1];
	}

	for(u64 i=0; i<bundles.size(); i++)
	{
		for(u64 j=0; j<bundles[i].point_idx.size(); j++)
		{
			// W block i, j
			Matrix A_block = A(bundles, points, bundles[i].point_idx[j], i);
			Matrix B_block = B(bundles, points, bundles[i].point_idx[j], i);

			Matrix W_block = transpose(A_block)*B_block;
	
			u64 stride = 45*(W_rows[i] + j);

			for(i64 c=0; c<3; c++)
				for(i64 l=0; l<15; l++)
					W_data[stride + c*15 + l] = W_block[l][c]; 
		}
	}

	W = SpMatrix((W_rows.size() - 1)*15, 3*bundles.size(), 15, 3, W_rows, W_cols_idx, W_data);
}

void hessian(std::vector<Bundle>& bundles, std::vector<Matrix>& points, SpMatrix& U, SpMatrix& V, SpMatrix& W)
{
	// Compute Hessian [[U, W], [W.T, V]]

	build_U(bundles, points, U);
	build_V(bundles, points, V);
	build_W(bundles, points, W);
}

}


