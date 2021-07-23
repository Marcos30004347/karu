#include <cmath>
#include <assert.h>
#include <iostream>

#include "camera/Camera.hpp"
#include "renderer/Renderer.hpp"

using namespace karu;

f32 radians(f32 degrees)
{
	return ( degrees * 3.14 ) / 180.f ;
}

int main()
{
	Renderer renderer;

	// 33mm focal lengths and (1000, 1000) principal point
	Camera cam = Camera(1187, 1187, 0, 0, algebra::Vector({0, 0, 0}), algebra::Vector({0, 0, radians(360)}), 10, 0.01, 0.0);

	std::vector<algebra::Vector> pixels;

	for(int i=-800; i<800; i+=20)
	{
		for(int j=-800; j<800; j+=20)
		{
			pixels.push_back(cam.projection({ (float)i, (float)j, 5000 }));
		}
	}
	
	for(algebra::Vector pix : pixels)
	{
		std::cout << (pix[0]) << " " << (pix[1]) << "\n";
	}

	renderer.draw2dPoints(pixels, 0, 0, 100, 100);

}
