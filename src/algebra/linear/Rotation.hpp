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
	// xy
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

	// 
	Matrix roll(3,3, {
		cos(a), -1*sin(a), 0,
		sin(a), cos(a), 0,
		0, 0, 1
	});

	return pitch*yaw*roll;
}

// Matrix rotationMaxtrixToAxisAngle(Matrix R)
// {
// 	Matrix u(3,1);

// 	f32 h = R[2][1];
// 	f32 f = R[1][2];
// 	f32 c = R[0][2];
// 	f32 g = R[2][0];
// 	f32 d = R[1][0];
// 	f32 b = R[0][1];

// 	u[0][0] = h - f;
// 	u[1][0] = c - g;
// 	u[2][0] = d - b;

// 	f32 theta = acos((R[0][0] + R[1][1] + R[2][2] - 1) / 2);

// 	u = u / (2 * sin(theta));

// 	std::cout << "u:\n";
// 	printMatrix(u/norm(u));

// 	return u*theta ;
// }

Matrix rotationMaxtrixToAxisAngle(Matrix m)
{
  f32 angle,x,y,z; // variables for result
	f32 epsilon = 0.01; // margin to allow for rounding errors
	f32 epsilon2 = 0.1; // margin to distinguish between 0 and 180 degrees

	// optional check that input is pure rotation, 'isRotationMatrix' is defined at:
	// https://www.euclideanspace.com/maths/algebra/matrix/orthogonal/rotation/
	if (
		(fabs(m[0][1]-m[1][0])< epsilon)
	  && (fabs(m[0][2]-m[2][0])< epsilon)
	  && (fabs(m[1][2]-m[2][1])< epsilon))
	{
		// singularity found
		// first check for identity matrix which must have +1 for all terms
		//  in leading diagonaland zero in other terms
		if (
			(fabs(m[0][1]+m[1][0]) < epsilon2)
		  && (fabs(m[0][2]+m[2][0]) < epsilon2)
		  && (fabs(m[1][2]+m[2][1]) < epsilon2)
		  && (fabs(m[0][0]+m[1][1]+m[2][2]-3) < epsilon2))
		{
			// this singularity is identity matrix so angle = 0
			return Matrix(3,1, {2*PI, 0, 0}); // zero angle, arbitrary axis
		}
		// otherwise this singularity is angle = 180
		angle = PI;
	
		f32 xx = (m[0][0]+1)/2;
		f32 yy = (m[1][1]+1)/2;
		f32 zz = (m[2][2]+1)/2;
		f32 xy = (m[0][1]+m[1][0])/4;
		f32 xz = (m[0][2]+m[2][0])/4;
		f32 yz = (m[1][2]+m[2][1])/4;
	
		if ((xx > yy) && (xx > zz)) 
		{ 
			// m[0][0] is the largest diagonal term
			if (xx < epsilon) 
			{
				x = 0;
				y = 0.7071;
				z = 0.7071;
			} else 
			{
				x = sqrt(xx);
				y = xy/x;
				z = xz/x;
			}
		} else if (yy > zz) { // m[1][1] is the largest diagonal term
			if (yy < epsilon) {
				x = 0.7071;
				y = 0;
				z = 0.7071;
			} else {
				y = sqrt(yy);
				x = xy/y;
				z = yz/y;
			}	
		} else { // m[2][2] is the largest diagonal term so base result on this
			if (zz < epsilon) {
				x = 0.7071;
				y = 0.7071;
				z = 0;
			} else {
				z = sqrt(zz);
				x = xz/z;
				y = yz/z;
			}
		}
	
		return Matrix(3,1, { x, y, z }) * angle; // return 180 deg rotation
	}

	// as we have reached here there are no singularities so we can handle normally
	f32 s = sqrt(
		 (m[2][1] - m[1][2])*(m[2][1] - m[1][2])
		+(m[0][2] - m[2][0])*(m[0][2] - m[2][0])
		+(m[1][0] - m[0][1])*(m[1][0] - m[0][1])
	); // used to normalise
	
	if (abs(s) < 0.001) s=1;

	// prevent divide by zero, should not happen if matrix is orthogonal and should be
	// caught by singularity test above, but I've left it in just in case
	angle = acos(( m[0][0] + m[1][1] + m[2][2] - 1)/2);

	x = (m[2][1] - m[1][2])/s;
	y = (m[0][2] - m[2][0])/s;
	z = (m[1][0] - m[0][1])/s;

	return Matrix(3,1, { x, y, z }) * angle;
}



Matrix axisAngleToRotationMaxtrix(Matrix a1) {

	f32 angle = norm(a1);
	a1 = a1/angle;

	double c = cos(angle);
	double s = sin(angle);
	double t = 1.0 - c;

	f32 m00 = c + a1[0][0]*a1[0][0]*t;
	f32 m11 = c + a1[1][0]*a1[1][0]*t;
	f32 m22 = c + a1[2][0]*a1[2][0]*t;

	f32 tmp1 = a1[0][0]*a1[1][0]*t;
	f32 tmp2 = a1[2][0]*s;
	f32 m10 = tmp1 + tmp2;
	f32 m01 = tmp1 - tmp2;

	tmp1 = a1[0][0]*a1[2][0]*t;
	tmp2 = a1[1][0]*s;

	f32 m20 = tmp1 - tmp2;
	f32 m02 = tmp1 + tmp2;    tmp1 = a1[1][0]*a1[2][0]*t;

	tmp2 = a1[0][0]*s;

	f32 m21 = tmp1 + tmp2;
	f32 m12 = tmp1 - tmp2;

	return Matrix(3,3, {
		m00, m01, m02,
		m10, m11, m12,
		m20, m21, m22,
	});
}


// Matrix axisAngleToRotationMaxtrix(Matrix u)
// {
// 	Matrix I(3,3);

// 	I[0][0] = 1;
// 	I[1][1] = 1;
// 	I[2][2] = 1;

// 	f32 angle = norm(u);

// 	Matrix axis = u / angle;

// 	Matrix skew(3,3, {
// 			0, -axis[2][0], axis[1][0],
// 			axis[2][0], 0, -axis[0][0],
// 			-axis[1][0], axis[0][0], 0,
// 	});

// 	Matrix outer_u(3,3, {
// 			axis[0][0]*axis[0][0], axis[1][0]*axis[0][0], axis[2][0]*axis[0][0],
// 			axis[0][0]*axis[1][0], axis[1][0]*axis[1][0], axis[2][0]*axis[1][0],
// 			axis[0][0]*axis[2][0], axis[1][0]*axis[2][0], axis[2][0]*axis[2][0],
// 	});

// 	return I*cos(angle) + skew*sin(angle) + outer_u*(1-cos(angle));
// }

} 
