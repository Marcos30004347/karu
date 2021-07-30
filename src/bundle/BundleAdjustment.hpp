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

#include <cmath>
#include <vector>
#include <set>
#include <array>
#include <assert.h>
#include "algebra/matrix/Matrix.hpp"
#include "algebra/sparse/SMatrix.hpp"
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

	for(Bundle b : bundles)
	{
		x.push_back(b.camera.fx);
		x.push_back(b.camera.fy);
		x.push_back(b.camera.cx);
		x.push_back(b.camera.cy);
		x.push_back(b.camera.Cx);
		x.push_back(b.camera.Cy);
		x.push_back(b.camera.Cz);
		x.push_back(b.camera.r1);
		x.push_back(b.camera.r2);
		x.push_back(b.camera.r3);
		x.push_back(b.camera.k1);
		x.push_back(b.camera.k2);
		x.push_back(b.camera.k3);
	}

	for(Point p : points)
	{
		x.push_back(p.x);
		x.push_back(p.y);
		x.push_back(p.z);
	}
}

// Point i, camera j
Matrix A(std::vector<Bundle>& bundles, std::vector<Point>& points, u64 i, u64 j)
{
	return Matrix(2, 13, {
		u_dfx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dfy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dcx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dcy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dCx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dCy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dCz(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dr1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dr2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dr3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dk1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dk2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dk3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dfx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dfy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dcx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dcy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dCx(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dCy(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dCz(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dr1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dr2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dr3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dk1(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dk2(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dk3(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
	});
}

Matrix B(std::vector<Bundle>& bundles, std::vector<Point>& points, u64 i, u64 j)
{
	return Matrix(2, 6, {
		u_dX(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dY(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		u_dZ(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dX(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dY(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
		v_dZ(bundles[j].camera.fx, bundles[j].camera.fy, bundles[j].camera.cx, bundles[j].camera.cy, bundles[j].camera.Cx, bundles[j].camera.Cy, bundles[j].camera.Cz, bundles[j].camera.r1, bundles[j].camera.r2, bundles[j].camera.r3, bundles[j].camera.k1, bundles[j].camera.k2, bundles[j].camera.k3, 0, 0, points[i].x, points[i].y, points[i].z),
	});
}

// Matrix denseHessian(std::vector<Bundle>& bundles, std::vector<Point>& points)
// {
// 	// Compute Hessian [[U, W], [W.T, V]]
	
// 	Matrix* U = new Matrix[bundles.size()]{ Matrix(2, 13) };
// 	Matrix* V = new Matrix[points.size()]{ Matrix(2, 6) };
// 	Matrix** W = new Matrix*[bundles.size()]{ nullptr };

// 	for(i64 j=0; j<bundles.size(); j++)
// 		W[j] = new Matrix[points.size()]{ Matrix(2, 6) };
	
// 	for(i64 j=0; j<bundles.size(); j++)
// 	{
// 		for(u64 i : bundles[j].point_idx)
// 		{
// 			U[j] = U[j] + A(bundles, points, i, j);
// 		}
// 	}

// 	for(i64 j=0; j<bundles.size(); j++)
// 	{
// 		for(u64 i : bundles[j].point_idx)
// 		{
// 			V[i] = V[i] + B(bundles, points, i, j);
// 		}
// 	}

// 	for(i64 j=0; j<bundles.size(); j++)
// 	{
// 		for(u64 i : bundles[j].point_idx)
// 		{
// 			W[i][j] = A(bundles, points, i, j) * B(bundles, points, i, j);
// 		}
// 	}
// }



void sparseHessian(std::vector<Bundle>& bundles, std::vector<Point>& points, SMatrix& U, SMatrix& V, SMatrix& W)
{
	// Compute Hessian [[U, W], [W.T, V]]
	std::vector<u64> U_rows;
	std::vector<u64> U_cols_idx;
	
	i64 U_col = 0;
	bool U_project_i = false;

	U_rows.push_back(0);
	U_cols_idx.push_back(0);

	for(i64 j=0; j<bundles.size(); j++)
	{
		U_rows.push_back(U_rows[j] + 1);
		U_cols_idx.push_back(U_cols_idx[j] + 13);
	}

	std::vector<f32> U_data(26*U_cols_idx.size(), 0);
	for(i64 j=0; j<bundles.size(); j++)
	{
	
		for(u64 i : bundles[j].point_idx)
		{
			Matrix U_block = A(bundles, points, i, j);
	
			// Push data into col major layout
			for(i64 c = 0; c<13; c++)
			{
				for(i64 l = 0; l<2; l++)
				{
					U_data[i*12 + l + c*2] += U_block[l][c];
				}
			}
		}
	}

	std::vector<u64> V_rows;
	std::vector<u64> V_cols_idx;

	V_rows.push_back(0);
	V_cols_idx.push_back(0);

	for(i64 i=0; i<points.size(); i++)
	{
		V_rows.push_back(V_rows[i] + 1);
		V_cols_idx.push_back(V_cols_idx[i] + 3);
	}

	std::vector<f32> V_data(6*V_cols_idx.size(), 0);
	for(i64 j=0; j<bundles.size(); j++)
	{
		for(u64 i : bundles[j].point_idx)
		{
			Matrix V_block = B(bundles, points, i, j);
	
			// Push data into col major layout
			for(i64 c = 0; c<3; c++)
			{
				for(i64 l = 0; l<2; l++)
				{
					V_data[i*6 + l + c*2] += V_block[l][c];
				}
			}
		}
	}
	
	U = SMatrix((U_rows.size() - 1)*2, 13*U_cols_idx.size(), 2, 13, U_rows, U_cols_idx, U_data);
	V = SMatrix((V_rows.size() - 1)*2, 6*V_cols_idx.size(), 2, 3, V_rows, V_cols_idx, V_data);
}

}


