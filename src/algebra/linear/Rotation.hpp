#pragma once

#include <cmath>
#include <tuple>

#include "algebra/matrix/Matrix.hpp"
#include "algebra/linear/Linear.hpp"

namespace karu::algebra {

Matrix lookAt(Matrix& from, Matrix& to)
{
	Matrix forward = from - to;
	forward = forward/norm(forward);

	Matrix tmp(3,1, { 0, 1, 0});

	Matrix right = cross(tmp, forward);
	right = right/norm(right);

	Matrix up = cross(forward, right);
	up = up/norm(up);

	Matrix lookAt(3,3);

	lookAt[0][0] = right[0][0];
	lookAt[0][1] = right[1][0];
	lookAt[0][2] = right[2][0];

	lookAt[1][0] = up[0][0];
	lookAt[1][1] = up[1][0];
	lookAt[1][2] = up[2][0];

	lookAt[2][0] = forward[0][0];
	lookAt[2][1] = forward[1][0];
	lookAt[2][2] = forward[2][0];

	return lookAt;
}

/**
 * Rotation over the x, y, z axes respectively
 */
Matrix getRotationMatrix(f32 y, f32 b, f32 a)
{

	Matrix pitch(3,3, {
		1, 0, 0,
		0, cos(y), -1*sin(y),
		0, sin(y), cos(y)
	});

	Matrix yaw(3,3, {
		cos(b), 0, sin(b),
		0, 1, 0,
		-1*sin(y), 0, cos(y)
	});

	Matrix roll(3,3, {
		cos(a), -1*sin(a), 0,
		sin(a), cos(a), 0,
		0, 0, 1
	});

	return pitch*yaw*roll;
}

Matrix rotationMaxtrixToAxisAngle(Matrix R)
{
	Matrix u(3,1);

	u[0][0] = R[2][1] - R[1][2];
	u[1][0] = R[0][2] - R[2][0];
	u[2][0] = R[1][0] - R[0][1];

	f32 angle = acos((R[0][0]+R[1][1]+R[2][2] - 1)/2);

	u = u/(2*sin(angle));

	return u*angle ;
}

Matrix axisAngleToRotationMaxtrix(Matrix u)
{
	Matrix I(3,3);

	I[0][0] = 1;
	I[1][1] = 1;
	I[2][2] = 1;

	f32 angle = norm(u);
	Matrix axis = u/angle;

	Matrix skew(3,3, {
			0, -1*axis[2][0], axis[1][0],
			axis[2][0], 0, -1*axis[0][0],
			-1*axis[1][0], axis[0][0], 0,
	});

	Matrix outer_u(3,3, {
			axis[0][0]*axis[0][0], axis[1][0]*axis[0][0], axis[2][0]*axis[0][0],
			axis[0][0]*axis[1][0], axis[1][0]*axis[1][0], axis[2][0]*axis[1][0],
			axis[0][0]*axis[2][0], axis[1][0]*axis[2][0], axis[2][0]*axis[2][0],
	});

	return I*cos(angle) + skew*sin(angle) + outer_u*(1-cos(angle));
}

} 
