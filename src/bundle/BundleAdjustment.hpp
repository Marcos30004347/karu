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
#include "algebra/linear/Linear.hpp"
#include "camera/CameraBundle.hpp"

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
	CameraBundle camera;

	// Camera observerd points pixel coordinates
	std::vector<Matrix> projections;

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
Matrix cameraParametersDerivatives(std::vector<Bundle>& bundles, std::vector<Matrix>& points, u64 i, u64 j)
{
	return Matrix(2, 15, {
		bund_u_dfx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dfy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dcx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dcy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dCx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dCy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dCz(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dr1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dr2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dr3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dk1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dk2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dk3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dp1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dp2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dfx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dfy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dcx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dcy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dCx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dCy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dCz(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dr1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dr2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dr3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dk1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dk2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dk3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dp1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dp2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
	});
}

Matrix pointParametersDerivatives(std::vector<Bundle>& bundles, std::vector<Matrix>& points, u64 i, u64 j)
{
	return Matrix(2, 3, {
		bund_u_dX(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dY(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_u_dZ(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),

		bund_v_dX(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dY(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
		bund_v_dZ(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.P[0][0], bundles[j].camera.P[1][0], bundles[j].camera.P[2][0], bundles[j].camera.R[0][0], bundles[j].camera.R[1][0], bundles[j].camera.R[2][0], bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i][0][0], points[i][1][0], points[i][2][0]),
	});
}

void buildHessianU(std::vector<Bundle>& bundles, std::vector<Matrix>& points, std::vector<std::vector<u64>>& point_idx_to_camera, f32 lambda, SpMatrix& U)
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
			Matrix A_block = cameraParametersDerivatives(bundles, points, i, j);
	
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
	
	// Augment the diagonal of U for Levenber-Marquardt
	if(lambda > 0.0)
	{
		for(i64 i=0; i<U.m_data.blocksCount(); i++)
		{
			for(i64 j=0; j<U.m_data.blockWidth(); j++)
			{
				u64 stride = i*U.m_data.blockHeight()*U.m_data.blockWidth() + j*U.m_data.blockHeight() + j;
				U.m_data.data()[stride] = U.m_data.data()[stride] + lambda;
			}
		}
	}
}

void buildHessianV(std::vector<Bundle>& bundles, std::vector<Matrix>& points, std::vector<std::vector<u64>>& point_idx_to_camera, f32 lambda, SpMatrix& V)
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
			Matrix B_block = pointParametersDerivatives(bundles, points, i, j);
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

	// Augment the diagonal of V for Levenber-Marquardt
	if(lambda > 0.0)
	{
		for(i64 i=0; i<V.m_data.blocksCount(); i++)
		{
			for(i64 j=0; j<V.m_data.blockWidth(); j++)
			{
				u64 stride = i*V.m_data.blockHeight()*V.m_data.blockWidth() + j*V.m_data.blockHeight() + j;
				V.m_data.data()[stride] = V.m_data.data()[stride] + lambda;
			}
		}
	}
}

void buildHessianVInverse(SpMatrix& V, SpMatrix& V_inv)
{
	std::vector<u64> V_rows = V.m_data.rowPtr();

	std::vector<u64> V_cols_idx = V.m_data.columnsIdx();

	std::vector<f32> V_data(V.m_data.storedElements(), 0);

	for(i64 i=0; i<V.m_data.blocksCount(); i++)
	{
		Matrix V_blk = Matrix(V.m_data.blockHeight(), V.m_data.blockWidth());
		
		for(i64 c=0; c<V.m_data.blockWidth(); c++)
			for(i64 l=0; l<V.m_data.blockHeight(); l++)
				V_blk[l][c] = V.m_data.data()[i*V.m_data.blockHeight()*V.m_data.blockWidth() + c*V.m_data.blockHeight() + l];

		std::pair<Matrix, Matrix> VLUP = LUPDecomposition(V_blk);
		Matrix V_block_inv = LUPInverse(VLUP.first, VLUP.second);

		for(i64 c=0; c<V.m_data.blockWidth(); c++)
			for(i64 l=0; l<V.m_data.blockHeight(); l++)
				V_data[i*V.m_data.blockHeight()*V.m_data.blockWidth() + c*V.m_data.blockHeight() + l] = V_block_inv[l][c];
	}

	V_inv = SpMatrix(V.m_data.lines(), V.m_data.columns(), V.m_data.blockHeight(), V.m_data.blockWidth(), V_rows, V_cols_idx, V_data);
}

void buildHessianUInverse(SpMatrix& U, SpMatrix& U_inv)
{
	std::vector<u64> U_rows = U.m_data.rowPtr();
	std::vector<u64> U_cols_idx = U.m_data.columnsIdx();
	std::vector<f32> U_data(U.m_data.storedElements(), 0);

	for(i64 i=0; i<U.m_data.blocksCount(); i++)
	{
		Matrix U_blk = Matrix(U.m_data.blockHeight(), U.m_data.blockWidth());
		
		for(i64 c=0; c<U.m_data.blockWidth(); c++)
			for(i64 l=0; l<U.m_data.blockHeight(); l++)
				U_blk[l][c] = U.m_data.data()[i*U.m_data.blockHeight()*U.m_data.blockWidth() + c*U.m_data.blockHeight() + l];

		std::pair<Matrix, Matrix> ULUP = LUPDecomposition(U_blk);
		Matrix U_block_inv = LUPInverse(ULUP.first, ULUP.second);

		for(i64 c=0; c<U.m_data.blockWidth(); c++)
			for(i64 l=0; l<U.m_data.blockHeight(); l++)
				U_data[i*U.m_data.blockHeight()*U.m_data.blockWidth() + c*U.m_data.blockHeight() + l] = U_block_inv[l][c];
	}

	U_inv = SpMatrix(U.m_data.lines(), U.m_data.columns(), U.m_data.blockHeight(), U.m_data.blockWidth(), U_rows, U_cols_idx, U_data);
}

void buildHessianW(std::vector<Bundle>& bundles, std::vector<Matrix>& points, std::vector<std::vector<u64>>& point_idx_to_camera, f32 lambda, SpMatrix& W)
{
	u64 W_block_cnt = 0;

	for(i64 i=0; i<bundles.size(); i++)
		for(u64 j : bundles[i].point_idx)
			W_block_cnt++;

	std::vector<f32> W_data(45*W_block_cnt, 0);
	std::vector<u64> W_cols_idx(W_block_cnt, 0);
	std::vector<u64> W_rows(bundles.size()+1, 0);

	// u64 idx = 0;

	// for(i64 i=0; i<points.size(); i++)
	// {
	// 	for(u64 j : point_idx_to_camera[i])
	// 	{
	// 		W_cols_idx[idx++] = 3*j;
	// 		W_rows[i+1]++;
	// 	}
	
	// 	if(i+2 < W_rows.size())
	// 		W_rows[i+2] = W_rows[i+1];
	// }

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

	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i=0; i<bundles[j].point_idx.size(); i++)
		{
			// W block i, j
			Matrix A_block = cameraParametersDerivatives(bundles, points, bundles[j].point_idx[i], j);
			Matrix B_block = pointParametersDerivatives(bundles, points, bundles[j].point_idx[i], j);

			Matrix W_block = transpose(A_block)*B_block;

			// block i j should be in j i block position
			u64 stride = 45*(W_rows[j] + i);

			for(i64 c=0; c<3; c++)
				for(i64 l=0; l<15; l++)
					W_data[stride + c*15 + l] = W_block[l][c]; 
		}
	}

	W = SpMatrix((W_rows.size() - 1)*15, 3*bundles.size(), 15, 3, W_rows, W_cols_idx, W_data);
}


void buildHessianW_T(std::vector<Bundle>& bundles, std::vector<Matrix>& points, std::vector<std::vector<u64>>& point_idx_to_camera, f32 lambda, SpMatrix& W_T)
{

	u64 W_block_cnt = 0;

	for(i64 i=0; i<bundles.size(); i++)
		for(u64 j : bundles[i].point_idx)
			W_block_cnt++;

	std::vector<f32> W_data(45*W_block_cnt, 0);
	std::vector<u64> W_cols_idx(W_block_cnt, 0);
	std::vector<u64> W_rows(points.size()+1, 0);

	u64 idx = 0;

	for(i64 i=0; i<points.size(); i++)
	{
		for(u64 j : point_idx_to_camera[i])
		{
			W_cols_idx[idx++] = 15*j;
			W_rows[i+1]++;
		}
	
		if(i+2 < W_rows.size())
			W_rows[i+2] = W_rows[i+1];
	}

	for(i64 i=0; i<points.size(); i++)
	{
		for(u64 j : point_idx_to_camera[i])
		{
			// W block i, j
			Matrix A_block = cameraParametersDerivatives(bundles, points, i, j);
			Matrix B_block = pointParametersDerivatives(bundles, points, i, j);

			Matrix W_block = transpose(A_block)*B_block;
			Matrix W = transpose(W_block);

			u64 stride = 45*(W_rows[i] + j);

			for(i64 c=0; c<15; c++)
				for(i64 l=0; l<3; l++)
					W_data[stride + c*3 + l] = W[l][c]; 
		}
	}

	W_T = SpMatrix((W_rows.size() - 1)*3, 15*bundles.size(), 3, 15, W_rows, W_cols_idx, W_data);
}


// Compute Hessian [[U, W], [W.T, V]]
void hessian(
	std::vector<Bundle>& bundles,
	std::vector<Matrix>& points,
	std::vector<std::vector<u64>>& point_idx_to_camera, f32 lambda,
	SpMatrix& U,
	SpMatrix& V,
	SpMatrix& W,
	SpMatrix& W_T
)
{
	buildHessianU(bundles, points, point_idx_to_camera, lambda, U);
	buildHessianV(bundles, points, point_idx_to_camera, lambda, V);
	buildHessianW(bundles, points, point_idx_to_camera, lambda, W);
	buildHessianW_T(bundles, points, point_idx_to_camera, lambda, W_T);
}

void shurComplementTimesX(SpMatrix& U, SpMatrix& W, SpMatrix& W_T, SpMatrix& V_inv, Matrix& x, Matrix& Sx)
{
	Matrix x1 = W_T*x;
	Matrix x2 = V_inv*x1;
	Matrix x3 = W*x2;
	Matrix x4 = U*x;
	Sx 				= x4 - x3;
}


void shurComplementConjugateGradient(SpMatrix& U, SpMatrix& W, SpMatrix& W_T, SpMatrix& V_inv, Matrix& x, Matrix b, float tolerance)
{
	Matrix Sx;
	Matrix Ap;

	shurComplementTimesX(U, W, W_T, V_inv, x, Sx);

	Matrix r = b - Sx;
	Matrix p = r;

	Matrix rsold = transpose(r)*r;

	u64 it = 0;

	while(norm(r) > tolerance && it++ < 20000)
	{
		// std::cout << "Converging: " << norm(r) << std::endl;
		
		shurComplementTimesX(U, W, W_T, V_inv, p, Ap);

		Matrix alpha = rsold/(transpose(p)*Ap);
	
		x = x + alpha*p;
		r = r - alpha*Ap;

		Matrix rsnew = transpose(r)*r;

		p = r + (rsnew/rsold)*p;
		rsold = rsnew;
	}
}


void solveNormalEquations(std::vector<Bundle>& bundles, std::vector<Matrix>& points, std::vector<std::vector<u64>>& point_idx_to_camera, f32& lambda)
{
	SpMatrix U, V, W, W_T, /*U_inv,*/ V_inv;

	std::vector<Matrix> r;
	std::vector<Matrix> rc;
	std::vector<Matrix> rp;
	
	lambda = 0;

	for(int j=0; j<bundles.size(); j++)
		r.push_back(Matrix(3, 1));
	
	for(int j=0; j<bundles.size(); j++)
	{
		rc.push_back(Matrix(15, 1));
	}
	
	for(int i=0; i<points.size(); i++)
		rp.push_back(Matrix(3, 1));


	for(u64 j=0; j<bundles.size(); j++)
	{

		for(u64 i=0; i<bundles[j].projections.size(); i++)
		{
			u64 k = bundles[j].point_idx[i];

			Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(points[k]);
			
			lambda += pow(norm(r),2);
			
			Matrix A_block = cameraParametersDerivatives(bundles, points, k, j);
			Matrix B_block = pointParametersDerivatives(bundles, points, k, j);

			rc[j] = rc[j] + transpose(A_block)*r;
			rp[k] = rp[k] + transpose(B_block)*r;
		}
	}

	lambda = sqrt(lambda);

	hessian(bundles, points, point_idx_to_camera, lambda, U, V, W, W_T);

	buildHessianVInverse(V, V_inv);
	// buildHessianUInverse(U, U_inv);

	for(i64 j=0; j<bundles.size(); j++)
	{
		Matrix Y_rp = Matrix(15, 1);
		for(u64 i : bundles[j].point_idx)
		{
			Matrix V_blk_inv(V_inv.m_data.blockHeight(), V_inv.m_data.blockWidth());
			
			for(i64 c=0; c<V_inv.m_data.blockWidth(); c++)
				for(i64 l=0; l<V_inv.m_data.blockHeight(); l++)
					V_blk_inv[l][c] = V_inv.m_data.data()[i*V_inv.m_data.blockHeight()*V_inv.m_data.blockWidth() + c*V_inv.m_data.blockHeight() + l];

			Matrix A_block = cameraParametersDerivatives(bundles, points, i, j);
			Matrix B_block = pointParametersDerivatives(bundles, points, i, j);

			Matrix Y_block = transpose(A_block)*B_block*V_blk_inv;

			Y_rp = Y_rp + Y_block*rp[i];
		}

		r[j] = rc[j] - Y_rp;
	}

	Matrix rf = Matrix(15*r.size(), 1);
	Matrix dC = Matrix(15*r.size(), 1);

	for(i64 k=0; k<r.size(); k++)
		for(i64 idx=0; idx<15; idx++)
			rf[k*15 + idx][0] = r[k][idx][0];

	shurComplementConjugateGradient(U, W, W_T, V_inv, dC, rf, 0.0000001f);

	Matrix dP = Matrix(3*points.size(), 1);

	for(i64 i=0; i<points.size(); i++)
	{
		Matrix K = Matrix(3, 1);
	
		for(i64 j: point_idx_to_camera[i])
		{
			Matrix A_block = cameraParametersDerivatives(bundles, points, bundles[j].point_idx[i], j);
			Matrix B_block = pointParametersDerivatives(bundles, points, bundles[j].point_idx[i], j);

			Matrix W_block = transpose(A_block)*B_block;
			Matrix W_block_T = transpose(W_block);
		
			Matrix dc(15, 1);
	
			for(i64 idx=0; idx<15; idx++)
			{
				dc[idx] = dC[j*15 + idx];
			}
	
			K = K + W_block_T*dc;
		}
	
		for(i64 idx=0; idx<3; idx++)
		{
			dP[i*3 + idx][0] = K[idx][0];
		}
	}

	// Update camera parameters
	for(i64 j=0; j<bundles.size(); j++)
	{
		bundles[j].camera.fx += dC[j*15 + 0][0];
		bundles[j].camera.fy += dC[j*15 + 1][0];
		bundles[j].camera.cx += dC[j*15 + 2][0];
		bundles[j].camera.cy += dC[j*15 + 3][0];
		bundles[j].camera.P[0][0] += dC[j*15 + 4][0];
		bundles[j].camera.P[1][0] += dC[j*15 + 5][0];
		bundles[j].camera.P[2][0] += dC[j*15 + 6][0];
		bundles[j].camera.R[0][0] += dC[j*15 + 7][0];
		bundles[j].camera.R[1][0] += dC[j*15 + 8][0];
		bundles[j].camera.R[2][0] += dC[j*15 + 9][0];
		bundles[j].camera.k1 += dC[j*15 + 10][0];
		bundles[j].camera.k2 += dC[j*15 + 11][0];
		bundles[j].camera.k3 += dC[j*15 + 12][0];
		bundles[j].camera.p1 += dC[j*15 + 13][0];
		bundles[j].camera.p2 += dC[j*15 + 14][0];
	}

	// Update point parameters
	for(i64 i=0; i<points.size(); i++)
	{
		points[i][0][0] += dP[i*3 + 0][0];
		points[i][1][0] += dP[i*3 + 1][0];
		points[i][2][0] += dP[i*3 + 2][0];
	}

	// for(i64 j=0; j<bundles.size(); j++)
	// {
	// 	std::cout << "fx: " << bundles[j].camera.fx << "\n";
	// 	std::cout << "fy: " << bundles[j].camera.fy << "\n";
	// 	std::cout << "cx: " << bundles[j].camera.cx << "\n";
	// 	std::cout << "cy: " << bundles[j].camera.cx << "\n";
	// 	std::cout << "Cx: " << bundles[j].camera.P[0][0] << "\n";
	// 	std::cout << "Cy: " << bundles[j].camera.P[1][0] << "\n";
	// 	std::cout << "Cz: " << bundles[j].camera.P[2][0] << "\n";
	// 	std::cout << "r1: " << bundles[j].camera.R[0][0] << "\n";
	// 	std::cout << "r2: " << bundles[j].camera.R[1][0] << "\n";
	// 	std::cout << "r3: " << bundles[j].camera.R[2][0] << "\n";
	// 	std::cout << "\n";
	// }
	// std::cout << "camera parameters:\n";
	// printMatrix(dC);
	// std::cout << "point parameters:\n";
	// printMatrix(dP);
}

}


