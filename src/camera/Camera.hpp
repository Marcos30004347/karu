#pragma once

#include "algebra/core/types.hpp"
#include "algebra/vector/Vector.hpp"
#include "camera/codegen/cameraModel.hpp"

namespace karu {


class Camera 
{
public:
	// focal lengths
	f32 fx, fy;

	// principal point
	f32 cx, cy;

	// Camera position
	f32 Cx, Cy, Cz;

	// angle and Rotation axis vector
	f32 r1, r2, r3;

	// distortion parameters
	f32 k1, k2, k3, p1, p2;
	Camera()
	{
		this->fx = 30;
		this->fy = 30;
		this->cx = 0;
		this->cy = 0;
		this->Cx = 0;
		this->Cy = 0;
		this->Cz = 0;
		this->r1 = 0;
		this->r2 = 0;
		this->r3 = ( 360 * 3.14 ) / 180.f;
		this->k1 = 0;
		this->k2 = 0;
		this->k3 = 0;
		this->p1 = 0;
		this->k2 = 0;
	}

	Camera(
		f32 fx, f32 fy,
		f32 cx, f32 cy, 
		algebra::Vector position_axis,
		algebra::Vector rotation_axis,
		f32 k1 = 0, f32 k2 = 0, f32 k3 = 0, f32 p1 = 0, f32 p2 = 0
	) 
	{
		f32 rx = rotation_axis[0];
		f32 ry = rotation_axis[1];
		f32 rz = rotation_axis[2];
		


		this->fx = fx;
		this->fy = fy;
		this->cx = cx;
		this->cy = cy;
		this->Cx = position_axis[0];
		this->Cy = position_axis[1];
		this->Cz = position_axis[2];
		this->r1 = rx;
		this->r2 = ry;
		this->r3 = rz;
		this->k1 = k1;
		this->k2 = k2;
		this->k3 = k3;
		this->p1 = p1;
		this->k2 = k2;
	}

	// Project the point (X,Y,Z)
	algebra::Vector projection(algebra::Vector pos)
	{
		f32 x = u(
			this->fx,
			this->fy,
			this->cx,
			this->cy,
			this->Cx,
			this->Cy,
			this->Cz,
			this->r1,
			this->r2,
			this->r3,
			this->k1,
			this->k2,
			this->k3,
			this->p1,
			this->p2,
			pos[0],
			pos[1],
			pos[2]
		);

		f32 y = v(
			this->fx,
			this->fy,
			this->cx,
			this->cy,
			this->Cx,
			this->Cy,
			this->Cz,
			this->r1,
			this->r2,
			this->r3,
			this->k1,
			this->k2,
			this->k3,
			this->p1,
			this->p2,
			pos[0],
			pos[1],
			pos[2]
		);
		return algebra::Vector({x, y});
	}


	algebra::Vector normalCoordinate(algebra::Vector pixel)
	{
		return { (pixel[0]-this->cx)/this->fx, (pixel[1]-this->cy)/this->fy };
	}

};

}
