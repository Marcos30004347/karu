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
#include "algebra/linear/Rotation.hpp"
#include "algebra/polynomial/Polynomial.hpp"
#include "camera/Camera.hpp"
#include "algebra/SVD/SVD.hpp"

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
	std::vector<Matrix> projections;

	// idx of the point that projections[i] is the projection in this camera
	std::vector<u64> point_idx;
};

f32 det3x3(Matrix& M)
{
	return 
		+ M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
		- M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
		+ M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
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

bool chooseRealizableSolution(Matrix Intrinsics[2], Matrix Rotations[2], Matrix Translations[2], Matrix* pts1, Matrix* pts2, i32 n, Matrix& Rotation, Matrix& Translation)
{
	i32 i, r, t, rotation = -1, translation = -1;

	Matrix cam0 = Intrinsics[0] * identity(3,4);
	Matrix cam1(3,4);
	
	bool realizable;

	for(r = 0; r < 2; r++)
	{
		realizable = true;

		for(t = 0; t < 2; t++)
		{
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
				rotation    = r;
				translation = t;
				break;
			}
		}
		if(realizable) break;
	}

	if(rotation != -1 && translation != -1)
	{
		Rotation    = Rotations[rotation];
		Translation = Translations[rotation];
	}

	return realizable;
}

void triangulatePoints(i32 n, Matrix* pts1, Matrix* pts2, Matrix cam0, Matrix cam1, std::vector<Matrix>& pts)
{
	for(i32 i = 0; i < n; i++)
	{
		pts.push_back(triangulate(pts1[i], pts2[i], cam0, cam1));
	}
}



Matrix skew(Matrix v)
{
	return Matrix(3,3, {
		0, -v[2][0], v[1][0],
		v[2][0], 0, -v[0][0],
		-v[1][0], v[0][0], 0,
	});
}

void estimateCameraFocalLengths(Matrix F, Matrix p1, Matrix p2, f32* f1, f32* f2)
{
	Matrix Uf, Vf_T;
	Matrix Uft, Vft_T;

	f32 Sf[3];
	f32 Sft[3];

	Matrix F_T = transpose(F); 

	svd(F, Uf, Sf, Vf_T);
	svd(F_T, Uft, Sft, Vft_T);

	Matrix e2  = Matrix(3, 1, { Vf_T[2][0], Vf_T[2][1], Vf_T[2][2] });

	Matrix e1 = Matrix(3, 1, { Vft_T[2][0], Vft_T[2][1], Vft_T[2][2] });

	Matrix I2 = identity(3,3);

	I2[2][2] = 0;

	Matrix f1_2_num = transpose(p1) * skew(e2) * I2 * F * p1 * transpose(p1) * transpose(F) * p1;
	Matrix f1_2_den = transpose(p1) * skew(e2) * I2 * transpose(F) * I2 * F * p1;
	
	Matrix f2_2_num = transpose(p2) * skew(e1) * I2 * transpose(F) * p2 * transpose(p2) * transpose(F) * p2;
	Matrix f2_2_den = transpose(p2) * skew(e1) * I2 * F * I2 * transpose(F) * p2;
	
	*f1 = f1_2_num[0][0]/f1_2_den[0][0];

	if(*f1 < 0)
	{
		*f1 = sqrt(-*f1);
	}
	else
	{
		*f1 = sqrt(*f1);
	}

	*f2 = f2_2_num[0][0]/f2_2_den[0][0];

	if(*f2 < 0)
	{
		*f2 = sqrt(-*f2);
	}
	else
	{
		*f2 = sqrt(*f2);
	}
}


void placeBundlesAndGetInitialPoints(
	std::vector<Bundle>& bundles,
	std::vector<Matrix>& points
)
{
	std::vector<std::vector<Matrix>> relativePosition;
	std::vector<std::vector<bool>>   relativeExists;

	for(i32 b = 0; b < bundles.size(); b++)
	{
		relativePosition.push_back(std::vector<Matrix>());
		relativeExists.push_back(std::vector<bool>());
	
		for(i32 k = 0; k < bundles.size(); k++)
			relativePosition[b].push_back(Matrix(3,4));
			relativeExists[b].push_back(false);
	}

	for(u32 b1 = 0; b1 < bundles.size(); b1++)
	{
		
		for(u32 b2 = 0; b2 < bundles.size(); b2++)
		{
		
			if(b1 == b2) continue;
		
			std::vector<std::pair<u64, u64>> matches;
		
			// Get matches
			for(u32 p1 = 0; p1 < bundles[b1].point_idx.size(); p1++)
			{
				for(u32 p2 = 0; p2 < bundles[b2].point_idx.size(); p2++)
				{
					if(bundles[b2].point_idx[p2] == bundles[b1].point_idx[p1])
					{
						matches.push_back({p1, p2});
					}
				}
			}
			
			if(matches.size() < 8)
			{
				// TODO: ignore bundle b2
				continue;
			}

			Matrix points1[8] = {
				bundles[b1].projections[matches[0].first],
				bundles[b1].projections[matches[1].first],
				bundles[b1].projections[matches[2].first],
				bundles[b1].projections[matches[3].first],
				bundles[b1].projections[matches[4].first],
				bundles[b1].projections[matches[5].first],
				bundles[b1].projections[matches[6].first],
				bundles[b1].projections[matches[7].first],
			};

			Matrix points2[8] = {
				bundles[b2].projections[matches[0].second],
				bundles[b2].projections[matches[1].second],
				bundles[b2].projections[matches[2].second],
				bundles[b2].projections[matches[3].second],
				bundles[b2].projections[matches[4].second],
				bundles[b2].projections[matches[5].second],
				bundles[b2].projections[matches[6].second],
				bundles[b2].projections[matches[7].second],
			};

			Matrix F = eightPointAlgorithm(points1, points2);

			f32 f1, f2;

			f32 a1 = (2*bundles[b1].camera.cx) / (2*bundles[b1].camera.cy);
			f32 a2 = (2*bundles[b2].camera.cx) / (2*bundles[b2].camera.cy);

			estimateCameraFocalLengths(
				F,
				Matrix(3,1, { bundles[b1].camera.cx, bundles[b1].camera.cy, 1 }),
				Matrix(3,1, { bundles[b2].camera.cx, bundles[b2].camera.cy, 1 }), 
				&f1,
				&f2
			);

			// Maybe need invertion of b1 by b2
			bundles[b1].camera.fx = a1 * f1;
			bundles[b1].camera.fy = f1;

			bundles[b2].camera.fx = a2 * f2;
			bundles[b2].camera.fy = f2;
			
			Matrix K1(3,3,
			{
				bundles[b1].camera.fx, 0, bundles[b1].camera.cx,
				0,    bundles[b1].camera.fy, bundles[b1].camera.cy,
				0,     0,   1
			});

			Matrix K2(3,3,
			{
				bundles[b2].camera.fx, 0, bundles[b2].camera.cx,
				0,    bundles[b2].camera.fy, bundles[b2].camera.cy,
				0,     0,   1
			});

			Matrix E = getEssentialMatrix(F, K1, K2);
			
			Matrix R1, R2, t1, t2;

			estimateRotationAndTranslation(E, R1, R2, t1, t2);

			Matrix Rs[2] = {R1, R2};
			Matrix Ts[2] = {t1, t2};
			Matrix Is[2] = {K1, K2};

			Matrix R, T;
		
			Matrix projections1[8] = {
				bundles[b1].projections[matches[0].first],
				bundles[b1].projections[matches[1].first],
				bundles[b1].projections[matches[2].first],
				bundles[b1].projections[matches[3].first],
				bundles[b1].projections[matches[4].first],
				bundles[b1].projections[matches[5].first],
				bundles[b1].projections[matches[6].first],
				bundles[b1].projections[matches[7].first],
			};

			Matrix projections2[8] = {
				bundles[b2].projections[matches[0].second],
				bundles[b2].projections[matches[1].second],
				bundles[b2].projections[matches[2].second],
				bundles[b2].projections[matches[3].second],
				bundles[b2].projections[matches[4].second],
				bundles[b2].projections[matches[5].second],
				bundles[b2].projections[matches[6].second],
				bundles[b2].projections[matches[7].second],
			};
	
			if(chooseRealizableSolution(Is, Rs, Ts, projections1, projections2, 8, R, T))
			{
				relativePosition[b1][b2] = Matrix(3, 4, {
					R[0][0], R[0][1], R[0][2], T[0][0],
					R[1][0], R[1][1], R[1][2], T[1][0],
					R[2][0], R[2][1], R[2][2], T[2][0],
				});

				relativeExists[b1][b2] = true;
			}
		}
	}

	for(i32 b1 = 0; b1 < bundles.size(); b1++)
	{
		for(i32 b2 = 0; b2 < bundles.size(); b2++)
		{
			if(b1 == b2) continue;

			if(relativeExists[b1][b2])
			{
				Matrix K1 = Matrix(3,3, {
					bundles[b1].camera.fx, 0, bundles[b1].camera.cx,
					0, bundles[b1].camera.fy, bundles[b1].camera.cy,
					0, 										 0,                     1,
				}) * identity(3,4);
		
				Matrix K2 = Matrix(3,3, {
					bundles[b1].camera.fx, 0, bundles[b1].camera.cx,
					0, bundles[b1].camera.fy, bundles[b1].camera.cy,
					0, 										 0,                     1,
				}) * relativePosition[b1][b2];
				
				triangulatePoints(8,  bundles[b1].projections.data(), bundles[b2].projections.data(), K1, K2, points);

				Matrix R(3,3, {
					relativePosition[b1][b2][0][0], relativePosition[b1][b2][0][1], relativePosition[b1][b2][0][2],
					relativePosition[b1][b2][1][0], relativePosition[b1][b2][1][1], relativePosition[b1][b2][1][2],
					relativePosition[b1][b2][2][0], relativePosition[b1][b2][2][1], relativePosition[b1][b2][2][2],
				});

				Matrix T(3,1, {
					relativePosition[b1][b2][0][3],
					relativePosition[b1][b2][1][3],
					relativePosition[b1][b2][2][3],
				});

				Camera camera1(
					bundles[b1].camera.fx,
					bundles[b1].camera.fy,
					bundles[b1].camera.cx,
					bundles[b1].camera.cy,
					Matrix(3,1, {0,0,0}),
					rotationMaxtrixToAxisAngle(getRotationMatrix(2*PI, 2*PI, 2*PI))
				);

				Camera camera2(
					bundles[b2].camera.fx,
					bundles[b2].camera.fy,
					bundles[b2].camera.cx,
					bundles[b2].camera.cy,
					transpose(R)*-1*T,
					rotationMaxtrixToAxisAngle(transpose(R))
				);

				bundles[b2].camera = camera2;
				bundles[b1].camera = camera1;
			}
		}
	}
}



// Point i, camera j
Matrix cameraParametersDerivatives(std::vector<Bundle>& bundles, std::vector<Matrix>& points, u64 i, u64 j)
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

Matrix pointParametersDerivatives(std::vector<Bundle>& bundles, std::vector<Matrix>& points, u64 i, u64 j)
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
			Matrix Aij = cameraParametersDerivatives(bundles, points, i, j);
	
			Matrix U_block = transpose(Aij) * Aij; 

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
	for(i64 i=0; i<U.m_data.blocksCount(); i++)
	{
		for(i64 j=0; j<U.m_data.blockWidth(); j++)
		{
			u64 idx = i*U.m_data.blockHeight()*U.m_data.blockWidth() + j*U.m_data.blockHeight() + j;
			U.m_data.data()[idx] = U.m_data.data()[idx] + lambda*U.m_data.data()[idx];
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
	for(i64 i=0; i<V.m_data.blocksCount(); i++)
	{
		for(i64 j=0; j<V.m_data.blockWidth(); j++)
		{
			u64 stride = i*V.m_data.blockHeight()*V.m_data.blockWidth() + j*V.m_data.blockHeight() + j;
			V.m_data.data()[stride] = V.m_data.data()[stride] + lambda*V.m_data.data()[stride];
		}
	}
}

void buildHessianVInverse(SpMatrix& V, SpMatrix& V_inv)
{
	std::vector<u64> V_rows 		= V.m_data.rowPtr();
	std::vector<f32> V_data			= std::vector<f32>(V.m_data.storedElements(), 0);
	std::vector<u64> V_cols_idx = V.m_data.columnsIdx();

	for(i64 i=0; i<V.m_data.blocksCount(); i++)
	{
		Matrix V_blk = Matrix(V.m_data.blockHeight(), V.m_data.blockWidth());
		
		for(i64 c=0; c<V.m_data.blockWidth(); c++)
		{
			for(i64 l=0; l<V.m_data.blockHeight(); l++)
			{
				V_blk[l][c] = V.m_data.data()[i*V.m_data.blockHeight()*V.m_data.blockWidth() + c*V.m_data.blockHeight() + l];
			}
		}

		std::pair<Matrix, Matrix> VLUP = LUPDecomposition(V_blk);
		Matrix V_block_inv             = LUPInverse(VLUP.first, VLUP.second);

		for(i64 c=0; c<V.m_data.blockWidth(); c++)
		{
			for(i64 l=0; l<V.m_data.blockHeight(); l++)
			{
				V_data[i*V.m_data.blockHeight()*V.m_data.blockWidth() + c * V.m_data.blockHeight() + l] = V_block_inv[l][c];
			}
		}
	}

	V_inv = SpMatrix(V.m_data.lines(), V.m_data.columns(), 3, 3, V_rows, V_cols_idx, V_data);
}

void buildHessianUInverse(SpMatrix& U, SpMatrix& U_inv)
{
	std::vector<u64> U_rows 		= U.m_data.rowPtr();
	std::vector<f32> U_data  		= std::vector<f32>(U.m_data.storedElements(), 0);
	std::vector<u64> U_cols_idx = U.m_data.columnsIdx();

	for(i64 i=0; i<U.m_data.blocksCount(); i++)
	{
		Matrix U_blk = Matrix(U.m_data.blockHeight(), U.m_data.blockWidth());
		
		for(i64 c=0; c<U.m_data.blockWidth(); c++)
		{
			for(i64 l=0; l<U.m_data.blockHeight(); l++)
			{
				U_blk[l][c] = U.m_data.data()[i*U.m_data.blockHeight()*U.m_data.blockWidth() + c*U.m_data.blockHeight() + l];
			}
		}

		std::pair<Matrix, Matrix> ULUP = LUPDecomposition(U_blk);
		Matrix U_block_inv 						 = LUPInverse(ULUP.first, ULUP.second);

		for(i64 c=0; c<U.m_data.blockWidth(); c++)
		{
			for(i64 l=0; l<U.m_data.blockHeight(); l++)
			{
				U_data[i*U.m_data.blockHeight()*U.m_data.blockWidth() + c*U.m_data.blockHeight() + l] = U_block_inv[l][c];
			}
		}
	}

	U_inv = SpMatrix(U.m_data.lines(), U.m_data.columns(), U.m_data.blockHeight(), U.m_data.blockWidth(), U_rows, U_cols_idx, U_data);
}

void buildHessianW(std::vector<Bundle>& bundles, std::vector<Matrix>& points, std::vector<std::vector<u64>>& point_idx_to_camera, f32 lambda, SpMatrix& W)
{
	u64 W_block_cnt = 0;

	for(i64 i=0; i < bundles.size(); i++)
	{
		W_block_cnt += bundles[i].point_idx.size();
	}

	std::vector<f32> W_data 		= std::vector<f32>(45*W_block_cnt, 0);
	std::vector<u64> W_cols_idx = std::vector<u64>(W_block_cnt, 0);
	std::vector<u64> W_rows 		= std::vector<u64>(bundles.size()+1, 0);

	u64 idx = 0;

	for(i64 i=0; i<bundles.size(); i++)
	{
		for(u64 j : bundles[i].point_idx)
		{
			W_rows[i+1] 			= W_rows[i+1] + 1;
			W_cols_idx[idx] = 3*j;
			idx = idx + 1;
		}

		if(i+2 < W_rows.size())
		{
			W_rows[i+2] = W_rows[i+1];
		}
	}

	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i : bundles[j].point_idx)
		{
			// W block i, j
			Matrix Aij = cameraParametersDerivatives(bundles, points, i, j);
			Matrix Bij = pointParametersDerivatives(bundles,  points, i, j);

			Matrix W_block = transpose(Aij)*Bij;

			for(i64 c=0; c<3; c++)
			{
				for(i64 l=0; l<15; l++)
				{
					W_data[45*(W_rows[j] + i) + c*15 + l] = W_block[l][c]; 
				}
			}
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

	while(norm(r) > tolerance && it < 2000)
	{
		shurComplementTimesX(U, W, W_T, V_inv, p, Ap);

		Matrix alpha = rsold/(transpose(p)*Ap);
	
		x = x + alpha*p;
		r = r - alpha*Ap;

		Matrix rsnew = transpose(r)*r;

		p = r + (rsnew/rsold)*p;
		rsold = rsnew;
		
		it++;
	}
}


void solveNormalEquations(std::vector<Bundle>& bundles, std::vector<Matrix>& points, std::vector<std::vector<u64>>& point_idx_to_camera, f32 tol, f32& lambda)
{
	SpMatrix U, V, W, W_T, /*U_inv,*/ V_inv;

	std::vector<Matrix> r;
	std::vector<Matrix> rc;
	std::vector<Matrix> rp;
	
	f32 error = 0;

	for(int j=0; j<bundles.size(); j++)
	{
		r.push_back(Matrix(3, 1));	
	}
	
	for(int j=0; j<bundles.size(); j++)
	{
		rc.push_back(Matrix(15, 1));
	}
	
	for(int i=0; i<points.size(); i++)
	{
		rp.push_back(Matrix(3, 1));
	}

	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i=0; i<bundles[j].projections.size(); i++)
		{
			u64 k = bundles[j].point_idx[i];

			Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(points[k]);

			error += pow(norm(r),2);

			Matrix Akj = cameraParametersDerivatives(bundles, points, k, j);
			Matrix Bkj = pointParametersDerivatives(bundles, points, k, j);

			rc[j] = rc[j] + transpose(Akj)*r;
			rp[k] = rp[k] + transpose(Bkj)*r;
		}
	}

	lambda = sqrt(error);

	hessian(bundles, points, point_idx_to_camera, lambda, U, V, W, W_T);

	buildHessianVInverse(V, V_inv);

	for(i64 j=0; j<bundles.size(); j++)
	{
		Matrix Y_rp = Matrix(15, 1);
		for(u64 i : bundles[j].point_idx)
		{
			Matrix V_blk_inv(V_inv.m_data.blockHeight(), V_inv.m_data.blockWidth());

			for(i64 c=0; c<V_inv.m_data.blockWidth(); c++)
			{
				for(i64 l=0; l<V_inv.m_data.blockHeight(); l++)
				{
					V_blk_inv[l][c] = V_inv.m_data.data()[i * V_inv.m_data.blockHeight() * V_inv.m_data.blockWidth() + c * V_inv.m_data.blockHeight() + l];
				}
			}

			Matrix A_block = cameraParametersDerivatives(bundles, points, i, j);
			Matrix B_block = pointParametersDerivatives(bundles, points, i, j);

			Matrix Y_block = transpose(A_block) * B_block * V_blk_inv;

			Y_rp = Y_rp + Y_block*rp[i];
		}

		r[j] = rc[j] - Y_rp;
	}

	Matrix rf = Matrix(15*r.size(), 1);
	Matrix dC = Matrix(15*r.size(), 1);

	for(i64 k=0; k<r.size(); k++)
	{
		for(i64 idx=0; idx<15; idx++)
		{
			rf[k*15 + idx][0] = r[k][idx][0];
		}
	}

	shurComplementConjugateGradient(U, W, W_T, V_inv, dC, rf, tol);

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
		
			K = K + W_block_T * Matrix(15,1, {
				dC[j*15 + 0][0],
				dC[j*15 + 1][0],
				dC[j*15 + 2][0],
				dC[j*15 + 3][0],
				dC[j*15 + 4][0],
				dC[j*15 + 5][0],
				dC[j*15 + 6][0],
				dC[j*15 + 7][0],
				dC[j*15 + 8][0],
				dC[j*15 + 9][0],
				dC[j*15 + 10][0],
				dC[j*15 + 11][0],
				dC[j*15 + 12][0],
				dC[j*15 + 13][0],
				dC[j*15 + 14][0],
			});
		}

		Matrix Vi_inv(V_inv.m_data.blockHeight(), V_inv.m_data.blockWidth());

		for(i64 c=0; c<V_inv.m_data.blockWidth(); c++)
		{
			for(i64 l=0; l<V_inv.m_data.blockHeight(); l++)
			{
				Vi_inv[l][c] = V_inv.m_data.data()[i * V_inv.m_data.blockHeight() * V_inv.m_data.blockWidth() + c * V_inv.m_data.blockHeight() + l];
			}
		}

		Matrix delta = Vi_inv*(rp[i] - K);

		dP[i*3 + 0][0] = delta[0][0];
		dP[i*3 + 1][0] = delta[1][0];
		dP[i*3 + 2][0] = delta[2][0];
	}

	// Update camera parameters
	// printMatrix(dC);
	for(i64 j=0; j<bundles.size(); j++)
	{
		// std::cout << dC[j*15 + 0][0] << " " << dC[j*15 + 1][0] << "\n";
		// std::cout << dC[j*15 + 2][0] << " " << dC[j*15 + 3][0] << "\n";

		bundles[j].camera.fx			+= dC[j*15 + 0][0];
		bundles[j].camera.fy 			+= dC[j*15 + 1][0];
		bundles[j].camera.cx 			+= dC[j*15 + 2][0];
		bundles[j].camera.cy 			+= dC[j*15 + 3][0];
		bundles[j].camera.P[0][0] += dC[j*15 + 4][0];
		bundles[j].camera.P[1][0] += dC[j*15 + 5][0];
		bundles[j].camera.P[2][0] += dC[j*15 + 6][0];
		bundles[j].camera.R[0][0] += dC[j*15 + 7][0];
		bundles[j].camera.R[1][0] += dC[j*15 + 8][0];
		bundles[j].camera.R[2][0] += dC[j*15 + 9][0];
		bundles[j].camera.k1 			+= dC[j*15 + 10][0];
		bundles[j].camera.k2 			+= dC[j*15 + 11][0];
		bundles[j].camera.k3 			+= dC[j*15 + 12][0];
		bundles[j].camera.p1 			+= dC[j*15 + 13][0];
		bundles[j].camera.p2 			+= dC[j*15 + 14][0];
	}

	// Update point parameters
	for(i64 i=0; i<points.size(); i++)
	{
		points[i][0][0] += dP[i*3 + 0][0];
		points[i][1][0] += dP[i*3 + 1][0];
		points[i][2][0] += dP[i*3 + 2][0];
		// std::cout << dP[i*3 + 0][0] << " "<< dP[i*3 + 1][0] << " " << dP[i*3 + 2][0] << "\n";
	}
}

void bundleAdjustment(std::vector<Bundle>& bundles, std::vector<Matrix>& wPoints, std::vector<std::vector<u64>> point_idx_to_bundle, f32 tol)
{
	placeBundlesAndGetInitialPoints(bundles, wPoints);

	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i=0; i<bundles[j].projections.size(); i++)
		{

			Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(wPoints[bundles[j].point_idx[i]]);
			std::cout << "Unadjusted Reprojection Error: " << norm(r) << "\n";
		}
	}

	f32 lambda = 1.0;
	i32 it = 0;
	i32 max_it = 20000;

	std::cout << "\n";

	do
	{
		it++;
		solveNormalEquations(bundles, wPoints, point_idx_to_bundle, lambda * lambda, lambda);
		std::cout << "Bundle Adjustment - step(" << it << "): error = " << lambda << std::endl;
	} while (lambda > tol && it < max_it);
	
	std::cout << "\n";
	std::cout << "Bundle adjustment takes " << it << " iterations\n";
	std::cout << "\n";

	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i=0; i<bundles[j].projections.size(); i++)
		{
				u64 k = bundles[j].point_idx[i];
				Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(wPoints[k]);
				std::cout << "Ajusted Reprojection Error: " << norm(r) << "\n";
		}
	}

	std::cout << "\n";

	for(i64 j=0; j<bundles.size(); j++)
	{
		std::cout << "Camera: " << j << "\n";

		std::cout << "fx: " << bundles[j].camera.fx << "\n";
		std::cout << "fy: " << bundles[j].camera.fy << "\n";
		std::cout << "cx: " << bundles[j].camera.cx << "\n";
		std::cout << "cy: " << bundles[j].camera.cx << "\n";
		std::cout << "Cx: " << bundles[j].camera.P[0][0] << "\n";
		std::cout << "Cy: " << bundles[j].camera.P[1][0] << "\n";
		std::cout << "Cz: " << bundles[j].camera.P[2][0] << "\n";
		std::cout << "r1: " << bundles[j].camera.R[0][0] << "\n";
		std::cout << "r2: " << bundles[j].camera.R[1][0] << "\n";
		std::cout << "r3: " << bundles[j].camera.R[2][0] << "\n";
		std::cout << "\n";
	}
}
}


