// #pragma once

// #include "bundle/BundleAdjustment.hpp"
// #include "camera/Camera.hpp"
// #include <random>
// #include <vector>
// #include <ctime>
// #include <tuple>

// using namespace karu;
// using namespace karu::algebra;


// algebra::Matrix randomVector()
// {
// 	std::mt19937 twister( std::time(0) ) ;
// 	std::uniform_real_distribution<f32> distr( -1000, 1000 ) ;

// 	algebra::Matrix vec(3,1, { distr(twister), distr(twister), distr(twister) });

// 	return algebra::Matrix(3,1, { vec[0]/norm(vec),  vec[1]/norm(vec),  vec[2]/norm(vec)});
// }

// std::vector<bundle::Bundle> genRandomBundles()
// {
// 	Camera cameras[25];

// 	algebra::Matrix points[25];
// 	algebra::Matrix dirs[25];
// 	algebra::Matrix positions[25];
// 	algebra::Matrix rotations[25];

// 	std::mt19937 dist_twister( std::time(0) ) ;
// 	std::uniform_real_distribution<f32> dist_distr( -1000, 1000 ) ;

// 	std::mt19937 point_twister( std::time(0) ) ;
// 	std::uniform_real_distribution<f32> point_distr( -100, 100 ) ;

// 	for(int i=0; i<25; i++)
// 		points[i] = algebra::Matrix({ point_distr(point_twister), point_distr(point_twister), point_distr(point_twister) });

// 	for(int i=0; i<25; i++)
// 	{
// 		dirs[i] = randomVector();
	
// 		f32 dist  = dist_distr(dist_twister);
// 		f32 angle = radians(360);
	
// 		positions[i] = algebra::Matrix({ dist*dirs[i][0], dist*dirs[i][1], dist*dirs[i][2] });
// 		rotations[i] = algebra::Matrix({  -1*dirs[i][0]*angle, -1*dirs[i][1]*angle,  -1*dirs[i][2]*angle });

// 		cameras[i] 	 = Camera(
// 			30, 30, 100, 100,
// 			positions[i],
// 			rotations[i],
// 			0,0,0,0,0
// 		);
// 	}

// 	std::vector<bundle::Bundle> bundles;

// 	for(int c=0;c<25; c++)
// 	{
	
// 		std::vector<u64> points_idx;
// 		std::vector<bundle::Pixel> projections;
	
// 		for(int i=0; i<25; i++)
// 		{
// 			algebra::Matrix pixel = cameras[c].projection(points[i]);
// 			algebra::Matrix norm_pixel = cameras[c].normalCoordinate(pixel);

// 			if(norm_pixel[0] < 1.0 && norm_pixel[1] < 1.0)
// 			{
// 				projections.push_back({ pixel[0][0], pixel[1][0] });
// 				points_idx.push_back(i);
// 			}
// 		}

// 		if(points_idx.size() > 4)
// 		{
// 			bundle::Bundle b;

// 			b.camera = cameras[c];
// 			b.point_idx = points_idx;
// 			b.projections = projections;

// 			bundles.push_back(b);
// 		}
// 	}	

// 	return bundles;
// }
