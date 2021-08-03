#include <cmath>
#include <assert.h>
#include <iostream>
#include <chrono>

#include "camera/Camera.hpp"
#include "renderer/Renderer.hpp"
#include "algebra/matrix/Matrix.hpp"
#include "algebra/linear/Rotation.hpp"

using namespace karu;

f32 radians(f32 degrees)
{
	return ( degrees * 3.14 ) / 180.f ;
}


std::chrono::steady_clock::time_point begin;
std::chrono::steady_clock::time_point end;

int main()
{
	f32 cy = 500;
	f32 cx = 500;

	f32 fy = 1000;
	f32 fx = 1000;

	Camera cam = Camera(fx, fy, cx, cy, algebra::Matrix(3, 1, {0, -2000, -2000}), algebra::Matrix(3, 1, { radians(45) , 0, 0 }), 0.5);

	std::vector<algebra::Matrix> pixels;

	for(int i=-1000; i<=1000; i+=20)
	{
		for(int j=-1000; j<=1000; j+=20)
		{
			begin = std::chrono::steady_clock::now();
		
			pixels.push_back(cam.projection(algebra::Matrix(3, 1, { (f32)i, (f32)j, 0 })));
		
			end = std::chrono::steady_clock::now();
		}
	}
	

	Renderer renderer(1000, 1000);
	renderer.draw2dPoints(cam, pixels);

}
