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

    Matrix tmp(3,1, {0,1,0});

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
}

Matrix getRotationMatrix(f32 yaw, f32 pitch, f32 roll)
{
    return Matrix(3, 3, {
        cos(yaw)*cos(pitch), cos(yaw)*sin(pitch)*sin(roll) - sin(yaw)*cos(roll), cos(yaw)*sin(pitch)*cos(roll) + sin(yaw)*sin(roll),
        sin(yaw)*cos(pitch), sin(yaw)*sin(pitch)*sin(roll) - cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll) - cos(yaw)*sin(roll),
        -1*sin(pitch),       cos(pitch)*sin(roll),                               cos(pitch)*cos(roll)
    });
}

std::tuple<Matrix, f32> rotationMaxtrixToAxisAngle(Matrix& R)
{
    Matrix u(3,1);

    u[0][0] = R[2][1] - R[1][2];
    u[1][0] = R[0][2] - R[2][0];
    u[2][0] = R[1][0] - R[0][1];
    
    u = u/norm(u);

    f32 angle = acos((R[0][0]+R[1][1]+R[2][2] - 1)/2);

    return { u, angle };
}

Matrix axisAngleToRotationMaxtrix(Matrix& axis, f32 angle)
{
    Matrix I(3,3);

    I[0][0] = 1;
    I[1][1] = 1;
    I[2][2] = 1;

    Matrix skew(3,3, {
        0, -1*axis[2][0], axis[1][0],
        axis[2][0], 0, -1*axis[0][0],
        -1*axis[1][0], axis[0][0], 0,
    });

    Matrix outer_u(3,3, {
        axis[0][0]*axis[0][0], axis[0][0]*axis[1][0], axis[0][0]*axis[2][0],
        axis[0][0]*axis[1][0], axis[1][0]*axis[1][0], axis[1][0]*axis[2][0],
        axis[0][0]*axis[2][0], axis[1][0]*axis[2][0], axis[2][0]*axis[2][0],
    });

    return cos(angle)*I + sin(angle)*skew + (1-cos(angle))*outer_u;
}

} 