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
	f32 cy = 3;
	f32 cx = 3;

	f32 fy = 30;
	f32 fx = 30;

	Camera cam = Camera(fx, fy, cx, cy, algebra::Vector({0, 0, 0}), algebra::Vector({0, 0, radians(360)}), 0.0);

	std::vector<algebra::Vector> pixels;

	for(int i=-400; i<400; i+=10)
	{
		for(int j=-400; j<400; j+=10)
		{
			begin = std::chrono::steady_clock::now();
			pixels.push_back(cam.projection({ (f32)i, (f32)j, 100 }));
			end = std::chrono::steady_clock::now();
			// std::cout << "Projection Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;
		}
	}
	
	for(algebra::Vector pix : pixels)
	{
		std::cout << (pix[0]-cx)/fx << " " << (pix[1]-cy)/fy << "\n";
	}

	renderer.draw2dPoints(pixels, cx, cy, fx, fy);

}
