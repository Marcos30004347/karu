#pragma once

#include "bundle/BundleAdjustment.hpp"
#include "camera/Camera.hpp"
#include <random>
#include <vector>
#include <ctime>
#include <tuple>

using namespace karu;
using namespace karu::algebra;


f32 norm(algebra::Vector v) {
	return std::sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] ) ;
}

algebra::Vector randomVector()
{
	std::mt19937 twister( std::time(0) ) ;
	std::uniform_real_distribution<f32> distr( -1000, 1000 ) ;

	algebra::Vector vec({ distr(twister), distr(twister), distr(twister) });

	return algebra::Vector({ vec[0]/norm(vec),  vec[1]/norm(vec),  vec[2]/norm(vec)});
}

f32 radians(f32 degrees)
{
	return ( degrees * 3.14 ) / 180.f ;
}


std::vector<bundle::Bundle> genRandomBundles()
{
	Camera cameras[25];

	algebra::Vector points[25];
	algebra::Vector dirs[25];
	algebra::Vector positions[25];
	algebra::Vector rotations[25];

	std::mt19937 dist_twister( std::time(0) ) ;
	std::uniform_real_distribution<f32> dist_distr( -1000, 1000 ) ;

	std::mt19937 point_twister( std::time(0) ) ;
	std::uniform_real_distribution<f32> point_distr( -100, 100 ) ;

	for(int i=0; i<25; i++)
		points[i] = algebra::Vector({ point_distr(point_twister), point_distr(point_twister), point_distr(point_twister) });

	for(int i=0; i<25; i++)
	{
		dirs[i] = randomVector();
	
		f32 dist  = dist_distr(dist_twister);
		f32 angle = radians(360);
	
		positions[i] = algebra::Vector({ dist*dirs[i][0], dist*dirs[i][1], dist*dirs[i][2] });
		rotations[i] = algebra::Vector({  -1*dirs[i][0]*angle, -1*dirs[i][1]*angle,  -1*dirs[i][2]*angle });

		cameras[i] 	 = Camera(
			30, 30, 100, 100,
			positions[i],
			rotations[i],
			0,0,0,0,0
		);
	}

	std::vector<bundle::Bundle> bundles;

	for(int c=0;c<25; c++)
	{
	
		std::vector<u64> points_idx;
		std::vector<bundle::Pixel> projections;
	
		for(int i=0; i<25; i++)
		{
			algebra::Vector pixel = cameras[c].projection(points[i]);
			algebra::Vector norm_pixel = cameras[c].normalCoordinate(pixel);

			if(norm_pixel[0] < 1.0 && norm_pixel[1] < 1.0)
			{
				projections.push_back({ pixel[0], pixel[1] });
				points_idx.push_back(i);
			}
		}

		if(points_idx.size() > 4)
		{
			bundle::Bundle b;

			b.camera = cameras[c];
			b.point_idx = points_idx;
			b.projections = projections;

			bundles.push_back(b);
		}
	}	

	return bundles;
}
