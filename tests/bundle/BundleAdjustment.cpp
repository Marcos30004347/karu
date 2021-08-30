#include <assert.h>

// #include "helpers.hpp"

#include "camera/Camera.hpp"
#include "renderer/Renderer.hpp"
#include "bundle/BundleAdjustment.hpp"
#include "algebra/linear/Rotation.hpp"
#include "algebra/compute/Compute.hpp"

using namespace karu;
using namespace karu::algebra;
using namespace karu::algebra::compute;
using namespace karu::bundle;


int main()
{
	// algebra::compute::Context::initContext();

	std::vector<Matrix> points;

	points.push_back(Matrix(3,1,{ 1.0,1.0,0.0 }));
	points.push_back(Matrix(3,1,{ 1.0,-1.0,0.0 }));
	points.push_back(Matrix(3,1,{ -1.0,-1.0,0.0 }));
	points.push_back(Matrix(3,1,{ -1.0,1.0,0.0 }));
	points.push_back(Matrix(3,1,{ 0.0, 0.0, sqrt(2.0) }));

	std::vector<Bundle> bundles;

	std::vector<Matrix> positions = {
		Matrix(3,1, { 0.0, -10.0, 1.0 }),
		Matrix(3,1, { 10.0, 0.0, 1.0 }),
		Matrix(3,1, { 0.0, 10.0, 1.0 }),
	};

	std::vector<Matrix> position_noises = {
		Matrix(3,1, { 8, 7, 5 }),
		Matrix(3,1, { 6, 4, 10 }),
		Matrix(3,1, { 9, 7, 1 }),
	};

	std::vector<Matrix> rotation_noises = {
		Matrix(3,1, { 0.1, 0.3, 0.1 }),
		Matrix(3,1, { 0.6, 0.7, 0.1 }),
		Matrix(3,1, { 0.3, 0.2, 0.1 }),
	};

	std::vector<Matrix> rotations = {
		Matrix(3,1, { 1.5707963 , 0.0, 0.0 }),
		Matrix(3,1, { 1.2091996, 1.2091996, 1.2091996 }),
		Matrix(3,1, { 0.0, 2.2214415, 2.2214415 }),
	};

	for(u64 i=0; i<3; i++)
	{
		Matrix R = axisAngleToRotationMaxtrix(rotations[i]);
		
		Matrix position = (transpose(R)*-1)*positions[i];
		Matrix rotation = rotationMaxtrixToAxisAngle(transpose(R));

		CameraBundle camera(
			3000,
			3000,
			500,
			500, 
			position,
			rotation
		);

		std::vector<u64>    point_idx;
		std::vector<Matrix> projections;

		for(int j=0; j<5; j++)
		{
			Matrix projection = camera.projection(points[j]);

			point_idx.push_back(j);
			projections.push_back(projection);
		}
	
		R = axisAngleToRotationMaxtrix(rotations[i]);
		
		position = (transpose(R)*-1)*positions[i] + position_noises[i];
		rotation = rotationMaxtrixToAxisAngle(transpose(R)) + rotation_noises[i];

		CameraBundle cam(
			3000,
			3000,
			500,
			500, 
			position,
			rotation
		);

		bundles.push_back({cam, projections, point_idx});
	}

	SpMatrix U, V, V_inv, W, W_T;

	std::vector<std::vector<u64>> point_idx_to_camera(points.size());

	for(i64 j=0; j<bundles.size(); j++)
		for(u64 i : bundles[j].point_idx)
			point_idx_to_camera[i].push_back(j);


	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i=0; i<bundles[j].projections.size(); i++)
		{
				u64 k = bundles[j].point_idx[i];
				Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(points[k]);
				std::cout << norm(r) << "\n";
		}
	}

	f32 lambda = 0.0;

	do
	{
		solveNormalEquations(bundles, points, point_idx_to_camera, lambda);
		std::cout << "converging: " << lambda << std::endl;
	} while (lambda > 0.1);
	
	
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "\n";

	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i=0; i<bundles[j].projections.size(); i++)
		{
				u64 k = bundles[j].point_idx[i];
				Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(points[k]);
				std::cout << "Reprojection Error: " << norm(r) << "\n";
		}
	}
	
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "\n";

	for(i64 j=0; j<bundles.size(); j++)
	{
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

	return 0;
}
