#include <assert.h>

#include "algebra/linear/Rotation.hpp"
#include "algebra/linear/Linear.hpp"

using namespace karu;
using namespace karu::algebra;

void should_rotate_points()
{
	Matrix R[3];

	R[0] = getRotationMatrix(0, 0, radians(90));
	R[1] = getRotationMatrix(0, radians(90), 0);
	R[2] = getRotationMatrix(radians(90), 0, 0);

	Matrix P0 = Matrix(3,1, { 1, 0, 0 });

	Matrix p0 = R[0]*P0;
	Matrix p1 = R[1]*P0;
	Matrix p2 = R[2]*P0;

	assert(roundToPrecision(p0[0][0], 3) == 0);
	assert(roundToPrecision(p0[0][1], 3) == 1);
	assert(roundToPrecision(p0[0][2], 3) == 0);

	assert(roundToPrecision(p1[0][0], 3) == 0);
	assert(roundToPrecision(p1[0][1], 3) == 0);
	assert(roundToPrecision(p1[0][2], 3) == -1);

	assert(roundToPrecision(p2[0][0], 3) == 1);
	assert(roundToPrecision(p2[0][1], 3) == 0);
	assert(roundToPrecision(p2[0][2], 3) == 0);
}

void should_convert_between_representations()
{
	Matrix R0 = Matrix(3,3, {
 		0.9983760, -0.0394619,  0.0410859,
   	0.0410859,  0.9983760, -0.0394619,
  	-0.0394619,  0.0410859,  0.9983760
	});

	Matrix axis_angle0 = rotationMaxtrixToAxisAngle(R0);
	
	assert(roundToPrecision(axis_angle0[0][0], 3) == 0.04);
	assert(roundToPrecision(axis_angle0[1][0], 3) == 0.04);
	assert(roundToPrecision(axis_angle0[2][0], 3) == 0.04);

	Matrix _R0 = axisAngleToRotationMaxtrix(axis_angle0);

	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			assert(roundToPrecision(_R0[i][j], 3) == roundToPrecision(R0[i][j], 3));
}

int main()
{
	should_rotate_points();
	should_convert_between_representations();

	// Look At
	// Matrix from = Matrix(3,1,{ 10, 10, 10});
	// Matrix to = Matrix(3,1,{0,0,0});

	// Matrix LookAt = lookAt(from, to);

	// for(int i=0; i<3; i++)
	// {
	// 	std::tuple<Matrix, f32> axis_angle = rotationMaxtrixToAxisAngle(LookAt);

	// 	Matrix axis = std::get<0>(axis_angle);
	// 	f32 angle = std::get<1>(axis_angle);
		
	// 	printMatrix(LookAt);
	// 	std::cout << "\n";
	// 	printMatrix(axis);
	// 	std::cout << "\n";

	// 	Matrix R_ = axisAngleToRotationMaxtrix(axis, angle);
	
	// 	printMatrix(R_);
	// 	std::cout << "\n";
	// }

	return 0;
}
