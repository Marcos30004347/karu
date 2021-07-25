#include <cmath>
#include <assert.h>
#include <iostream>
#include <chrono>

#include "camera/Camera.hpp"
#include "renderer/Renderer.hpp"

using namespace karu;

f32 radians(f32 degrees)
{
	return ( degrees * 3.14 ) / 180.f ;
}


std::chrono::steady_clock::time_point begin;
std::chrono::steady_clock::time_point end;

int main()
{
	Renderer renderer(1000, 1000);

	// 33mm focal lengths and (1000, 1000) principal point
	Camera cam = Camera(1187, 1187, 100, 100, algebra::Vector({0, 0, 0}), algebra::Vector({0, 0, radians(360)}), 15.3);

	std::vector<algebra::Vector> pixels;

	for(int i=-800; i<800; i+=10)
	{
		for(int j=-800; j<800; j+=10)
		{
			begin = std::chrono::steady_clock::now();
			pixels.push_back(cam.projection({ (float)i, (float)j, 5000 }));
			end = std::chrono::steady_clock::now();
			std::cout << "Projection Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;
		}
	}
	
	// for(algebra::Vector pix : pixels)
	// {
	// 	std::cout << (pix[0]) << " " << (pix[1]) << "\n";
	// }

	renderer.draw2dPoints(pixels, 100, 100, 100, 100);

}
