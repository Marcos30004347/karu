#include <assert.h>

#include "camera/Camera.hpp"
#include "camera/CameraBundle.hpp"
#include "algebra/linear/Rotation.hpp"

using namespace karu;
using namespace karu::algebra;


int main()
{
	Matrix p1 = Matrix(3,1,{ 1.0, 1.0,0.0 });

	std::vector<Matrix> positions = {
		Matrix(3,1, { 0.0, -10.0, 1.0 }),
		Matrix(3,1, { 10.0, 0.0, 1.0 }),
		Matrix(3,1, { 0.0, 10.0, 1.0 }),
	};

	std::vector<Matrix> rotations = {
		Matrix(3,1, { 1.5707963 , 0.0, 0.0 }),
		Matrix(3,1, { 1.2091996, 1.2091996, 1.2091996 }),
		Matrix(3,1, { 0.0, 2.2214415, 2.2214415 }),
	};

	Camera camera(
		3000,
		3000,
		500,
		500, 
		positions[2],
		rotations[2]
	);

	Matrix R = axisAngleToRotationMaxtrix(rotations[2]);

	CameraBundle camera_bundle(
		3000,
		3000,
		500,
		500, 
		(transpose(R)*-1)*positions[2],
		rotationMaxtrixToAxisAngle(transpose(R))
	);

	Matrix p1_proj_cb = camera_bundle.projection(p1);
	Matrix p1_proj_c  = camera.projection(p1);

	assert(p1_proj_cb[0][0] == p1_proj_c[0][0]);

	return 0;
}
