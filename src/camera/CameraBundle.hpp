#pragma once

#include "algebra/core/types.hpp"
#include "algebra/vector/Vector.hpp"
#include "algebra/matrix/Matrix.hpp"
#include "algebra/linear/Linear.hpp"
#include "camera/codegen/CameraModelBundle.hpp"

namespace karu {

class CameraBundle 
{
public:
	// focal lengths
	f32 fx, fy;

	// principal point
	f32 cx, cy;

	// Camera position
	algebra::Matrix P;

	// angle and Rotation axis vector
	algebra::Matrix R;

	// distortion parameters
	f32 k1, k2, k3, p1, p2;

	CameraBundle()
	{
		this->R = algebra::Matrix(3,1, {0,0,2*PI});
		this->P = algebra::Matrix(3,1, {0,0,0});
	
		this->fx = -1;
		this->fy = -1;
		this->cx = -1;
		this->cy = -1;
	
		this->k1 = 0;
		this->k2 = 0;
		this->k3 = 0;
		this->p1 = 0;
		this->k2 = 0;
	}

	CameraBundle(
		f32 fx, f32 fy,
		f32 cx, f32 cy, 
		algebra::Matrix P,
		algebra::Matrix R,
		f32 k1 = 0, f32 k2 = 0, f32 k3 = 0,
		f32 p1 = 0, f32 p2 = 0
	) 
	{
		this->R = R;
		this->P = P;

		this->fx = fx;
		this->fy = fy;
		this->cx = cx;
		this->cy = cy;
		this->k1 = k1;
		this->k2 = k2;
		this->k3 = k3;
		this->p1 = p1;
		this->k2 = k2;
	}

	// Project the point (X,Y,Z)
	algebra::Matrix projection(algebra::Matrix pos)
	{
		f32 x = bund_u(
			this->fx,
			this->fy,
			this->cx,
			this->cy,
			this->P[0][0],
			this->P[1][0],
			this->P[2][0],
			this->R[0][0],
			this->R[1][0],
			this->R[2][0],
			this->k1,
			this->k2,
			this->k3,
			this->p1,
			this->p2,
			pos[0][0],
			pos[1][0],
			pos[2][0]
		);

		f32 y = bund_v(
			this->fx,
			this->fy,
			this->cx,
			this->cy,
			this->P[0][0],
			this->P[1][0],
			this->P[2][0],
			this->R[0][0],
			this->R[1][0],
			this->R[2][0],
			this->k1,
			this->k2,
			this->k3,
			this->p1,
			this->p2,
			pos[0][0],
			pos[1][0],
			pos[2][0]
		);
	
		return algebra::Matrix(3,1, {x, y, 1});
	}
};

}
