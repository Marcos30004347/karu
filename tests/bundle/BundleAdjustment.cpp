#include <assert.h>

// #include "helpers.hpp"

#include "camera/Camera.hpp"
#include "renderer/Renderer.hpp"
#include "bundle/BundleAdjustment.hpp"
#include "algebra/linear/Rotation.hpp"

using namespace karu;
using namespace karu::algebra;
using namespace karu::bundle;


int main()
{
	std::vector<Matrix> points;

	points.push_back(Matrix(3,1,{ 1.0,1.0,0.0 }));
	points.push_back(Matrix(3,1,{ 1.0,-1.0,0.0 }));
	points.push_back(Matrix(3,1,{ -1.0,-1.0,0.0 }));
	points.push_back(Matrix(3,1,{ -1.0,1.0,0.0 }));
	points.push_back(Matrix(3,1,{ 0.0, 0.0, sqrt(2.0) }));

	std::vector<Bundle> bundles(3);

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

	Matrix O = Matrix(3,1, {0,0,0});

	for(u64 i=0; i<3; i++)
	{
		bundles[i].camera = Camera(
			3000,
			3000,
			500,
			500, 
			positions[i],
			rotations[i]
		);
		
		std::vector<Matrix> pixels;

		for(int j=0; j<5; j++)
		{
			Matrix projection = bundles[i].camera.projection(points[j]);
		
			bundle::Pixel p;
			p.u = projection[0][0];
			p.v = projection[1][0];
		
			bundles[i].point_idx.push_back(j);
			bundles[i].projections.push_back(p);

			pixels.push_back(projection);
		}

		// Renderer renderer(1000, 1000);
		// renderer.draw2dPoints(bundles[i].camera, pixels);
	}

	SpMatrix U, V, W;

	hessian(bundles, points, U, V, W);
	printMatrix(U);
	std::cout << "\n";
	printMatrix(V);
	std::cout << "\n";
	printMatrix(W);
	return 0;
}
