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

int shouldAjustAlmostCompletedBundles()
{
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
		Matrix(3,1, { 0.005, 0.04, 0.01 }),
		Matrix(3,1, { 0.01, 0.006, 0.05 }),
		Matrix(3,1, { 0.004, 0.01, 0.003 }),
	};

	std::vector<Matrix> rotation_noises = {
		Matrix(3,1, { 0.001, 0.03, 0.01 }),
		Matrix(3,1, { 0.06, 0.007, 0.01 }),
		Matrix(3,1, { 0.03, 0.02, 0.01 }),
	};

	std::vector<Matrix> rotations = {
		Matrix(3,1, { 1.5707963 , 0.0, 0.0 }),
		Matrix(3,1, { 1.2091996, 1.2091996, 1.2091996 }),
		Matrix(3,1, { 0.0, 2.2214415, 2.2214415 }),
	};

	for(u64 i=0; i<3; i++)
	{
		Camera camera(
			1520,
			1520,
			500,
			500, 
			positions[i],
			rotations[i]
		);

		std::vector<u64>    point_idx;
		std::vector<Matrix> projections;

		for(int j=0; j<5; j++)
		{
			Matrix projection = camera.projection(points[j]);
			point_idx.push_back(j);
			projections.push_back(projection);
		}

		Camera cam(
			1520,
			1520,
			500,
			500, 
			positions[i] + position_noises[i],
			rotations[i] + rotation_noises[i]
		);

		bundles.push_back({cam, projections, point_idx});
	}

	SpMatrix U, V, V_inv, W, W_T;

	std::vector<std::vector<u64>> point_idx_to_camera(points.size());

	for(i64 j=0; j<bundles.size(); j++)
		for(u64 i : bundles[j].point_idx)
			point_idx_to_camera[i].push_back(j);

	points[0] = points[0] + Matrix(3,1, {0.1, 0.01, 0.01});
	points[1] = points[1] + Matrix(3,1, {0.01, 0.01, 0.01});
	points[2] = points[2] + Matrix(3,1, {0.01, 0.0, 0.01});
	points[3] = points[3] + Matrix(3,1, {0.01, 0.01, 0.1});
	points[4] = points[4] + Matrix(3,1, {0.01, 0.1, 0.01});

	// for(u64 j=0; j<bundles.size(); j++)
	// {
	// 	for(u64 i=0; i<bundles[j].projections.size(); i++)
	// 	{
	// 			Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(points[bundles[j].point_idx[i]]);
	// 			std::cout << norm(r) << "\n";
	// 	}
	// }

	f32 lambda = 1.0;
	i32 it = 0;

	do
	{
		it++;
		solveNormalEquations(bundles, points, point_idx_to_camera, lambda * lambda, lambda);
		// std::cout << "step(" << it << "): error = " << lambda << std::endl;
	} while (lambda > 0.9 && it < 2000);
	
	std::cout << "\n";

	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i=0; i<bundles[j].projections.size(); i++)
		{
				u64 k = bundles[j].point_idx[i];
				Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(points[k]);
				assert(norm(r) < 0.9);
				// std::cout << "Reprojection Error: " << norm(r) << "\n";
		}
	}

	// std::cout << "\n";

	// for(i64 j=0; j<bundles.size(); j++)
	// {
	// 	std::cout << "fx: " << bundles[j].camera.fx << "\n";
	// 	std::cout << "fy: " << bundles[j].camera.fy << "\n";
	// 	std::cout << "cx: " << bundles[j].camera.cx << "\n";
	// 	std::cout << "cy: " << bundles[j].camera.cx << "\n";
	// 	std::cout << "Cx: " << bundles[j].camera.P[0][0] << "\n";
	// 	std::cout << "Cy: " << bundles[j].camera.P[1][0] << "\n";
	// 	std::cout << "Cz: " << bundles[j].camera.P[2][0] << "\n";
	// 	std::cout << "r1: " << bundles[j].camera.R[0][0] << "\n";
	// 	std::cout << "r2: " << bundles[j].camera.R[1][0] << "\n";
	// 	std::cout << "r3: " << bundles[j].camera.R[2][0] << "\n";
	// 	std::cout << "\n";
	// }

	return 0;
}



int shouldBundleAdjustIncompletedBundles()
{
	std::vector<Matrix> points;

	points.push_back(Matrix(3,1,{ -1.0, -1.5, 0.0 }));
	points.push_back(Matrix(3,1,{ 1.7, 1., 0.0 }));
	points.push_back(Matrix(3,1,{ .7, .5, 0.0 }));
	points.push_back(Matrix(3,1,{ -0.5, -1.1, 0.0 }));
	points.push_back(Matrix(3,1,{ -1.0, -1.5, .5 }));
	points.push_back(Matrix(3,1,{ 1.7, 1., 1.0 }));
	points.push_back(Matrix(3,1,{ .7, .5, 1.0 }));
	points.push_back(Matrix(3,1,{ -0.5, -1.1, 1.0 }));

	std::vector<Matrix> positions = {
		Matrix(3,1, { 0.0, -15.0, 1.0 }),
		Matrix(3,1, { 15.0, 0.0, 1.0 }),
		Matrix(3,1, { 15.0, -15.0, 1.0 }),
		Matrix(3,1, { 0.0, 0.0, 15.0 }),
	};

	std::vector<Matrix> rotations = {
		Matrix(3,1, { 1.5707963 , 0.0, 0.0 }),
		Matrix(3,1, { 1.2091996, 1.2091996, 1.2091996 }),
		Matrix(3,1, { 1.4821898, 0.6139431, 0.6139431 }),
		Matrix(3,1, { 0, 0, 1.5707963 }),
	};

	std::vector<Matrix> projections[2];

	for(u64 i=0; i<2; i++)
	{
		Camera camera(
			1600,
			1600,
			500,
			500,
			positions[i],
			rotations[i]
		);

		std::vector<u64> point_idx;

		for(int j=0; j<points.size(); j++)
		{
			Matrix projection = camera.projection(points[j]);
			projections[i].push_back(projection);
			point_idx.push_back(j);
		}
	}

	std::vector<Bundle> bundles;
	std::vector<Matrix> wPoints;

	bundles.push_back({
		Camera(
			0,
			0,
			500,
			500,
			Matrix(3,1, {0, 0, 0}),
			rotationMaxtrixToAxisAngle(getRotationMatrix(2*PI, 2*PI, 2*PI))
		),
		{	
			projections[0][0],
			projections[0][1],
			projections[0][2],
			projections[0][3],
			projections[0][4],
			projections[0][5],
			projections[0][6],
			projections[0][7]
		},
		{0, 1, 2, 3, 4, 5, 6, 7}
	});

	bundles.push_back({
		Camera(
			0,
			0,
			500,
			500,
			Matrix(3,1, {0, 0, 0}),
			rotationMaxtrixToAxisAngle(getRotationMatrix(2*PI, 2*PI, 2*PI))
		),
		{	
			projections[1][0],
			projections[1][1],
			projections[1][2],
			projections[1][3],
			projections[1][4],
			projections[1][5],
			projections[1][6],
			projections[1][7]
		},
		{0, 1, 2, 3, 4, 5, 6, 7}
	});

	std::vector<std::vector<u64>> point_idx_to_camera(points.size());

	for(i64 j=0; j<bundles.size(); j++)
		for(u64 i : bundles[j].point_idx)
			point_idx_to_camera[i].push_back(j);

	bundleAdjustment(bundles, wPoints, point_idx_to_camera, 0.000001);

	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i=0; i<bundles[j].projections.size(); i++)
		{

			Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(wPoints[bundles[j].point_idx[i]]);
			assert(norm(r) < 0.000001);
		}
	}

	return 0;
}


int main()
{
	// shouldAjustAlmostCompletedBundles();
	shouldBundleAdjustIncompletedBundles();
}
