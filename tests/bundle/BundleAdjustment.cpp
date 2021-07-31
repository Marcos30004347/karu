#include <assert.h>

#include "helpers.hpp"

#include "bundle/BundleAdjustment.hpp"
#include "camera/Camera.hpp"
#include "renderer/Renderer.hpp"

using namespace karu;
using namespace karu::algebra;
using namespace karu::bundle;

Vector cross(Vector A, Vector B)
{
	return {
		A[1]*B[2] - A[2]*B[1],
		A[2]*B[0] - A[0]*B[2],
		A[0]*B[1] - A[1]*B[0],
	};
}
 
std::vector<Bundle> genBundle()
{
	Vector points[4] = {
		{1.0,0.0,0.0},
		{0.0,1.0,0.0},
		{0.0,0.0,1.0},
		{0.0,0.0,0.0},
		// {0.5,0.5,0.5},
		// {0.5,-0.5,0.5},
		// {-0.5,0.5,0.5},
		// {1.0,1.0,0.5},
	};

	Vector positions[] = {
		{10, 10, 10},
		{10, 5, 5},
		{5, 10, 5},
	};

	Vector rays[] = {
		cross({-10, -10, -10}, {0, 0, 1}),
		{-2, -1, -1},
		{-1, -2, -1},
	};

	std::cout << rays[0][0] << " " << rays[0][1] << " " << rays[0][2] << "\n";
	rays[0] = Vector({rays[0][0]/norm(rays[0]), rays[0][1]/norm(rays[0]), rays[0][2]/norm(rays[0])});
	// rays[1] = Vector({rays[1][0]/norm(rays[1]), rays[1][1]/norm(rays[1]), rays[1][2]/norm(rays[1])});
	// rays[2] = Vector({rays[2][0]/norm(rays[2]), rays[2][1]/norm(rays[2]), rays[2][2]/norm(rays[2])});

	std::cout << rays[0][0] << " " << rays[0][1] << " " << rays[0][2] << "\n";
	
	rays[0] = Vector({rays[0][0]*radians(45), rays[0][1]*radians(45), rays[0][2]*radians(45)});
	// rays[1] = Vector({rays[1][0]*radians(360), rays[1][1]*radians(360), rays[1][2]*radians(360)});
	// rays[2] = Vector({rays[2][0]*radians(360), rays[2][1]*radians(360), rays[2][2]*radians(360)});

	std::cout << rays[0][0] << " " << rays[0][1] << " " << rays[0][2] << "\n";

	Bundle bundles[3];

	bundles[0].camera = Camera(1, 1, 0, 0, positions[0], rays[0]);
	bundles[1].camera = Camera(1, 1, 0, 0, positions[1], rays[1]);
	bundles[2].camera = Camera(1, 1, 0, 0, positions[2], rays[2]);


	std::vector<Vector> pixels;
	for(int i=0; i<4; i++)
	{
		Vector projection = bundles[0].camera.projection(points[i]);
		
		pixels.push_back(projection);
	
		std::cout
			<< (projection[0]/bundles[0].camera.fx)
			<< " "
			<< (projection[1]/bundles[0].camera.fx)
			<< "\n";
	}


	Renderer renderer(1000, 1000);

	renderer.draw2dPoints(pixels, 0, 0, bundles[0].camera.fx, bundles[0].camera.fy);

	return { bundles[0], bundles[1], bundles[2] };
}

int main()
{
	genBundle();
	// std::vector<Bundle> bundles = genRandomBundles();

	// std::cout << bundles.size() << std::endl;

	return 0;
}
