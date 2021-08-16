#include <assert.h>

#include "bundle/Homography.hpp"
#include "bundle/BundleAdjustment.hpp"

using namespace karu;
using namespace karu::bundle;
using namespace karu::algebra;


int main()
{
	std::vector<Matrix> points;

	points.push_back(Matrix(3,1,{ 1.0, 1.0, 0.0 }));
	points.push_back(Matrix(3,1,{ -1.0, 1.0, 0.0 }));
	points.push_back(Matrix(3,1,{ -1.0, -1.0, 0.0 }));
	points.push_back(Matrix(3,1,{ 1.0, -1.0, 0.0 }));
	points.push_back(Matrix(3,1,{ 1.0, 1.0, 1.0 }));
	points.push_back(Matrix(3,1,{ -1.0, 1.0, 1.0 }));
	points.push_back(Matrix(3,1,{ -1.0, -1.0, 1.0 }));
	points.push_back(Matrix(3,1,{ 1.0, -1.0, 1.0 }));

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
		Camera camera(
			3000,
			3000,
			500,
			500, 
			positions[i],
			rotations[i]
		);
	
		std::vector<u64> point_idx;
		std::vector<Matrix> projections;

		for(int j=0; j<points.size(); j++)
		{
			Matrix projection = camera.projection(points[j]);
			point_idx.push_back(j);
			projections.push_back(projection);
		}

		Camera cam(
			3000,
			3000,
			500,
			500, 
			positions[i],
			rotations[i]
		);

		bundles.push_back({cam, projections, point_idx});
	}

	// normalize points
	for(u64 i=0; i<3; i++)
	{
		Matrix m(3,1);

		
		for(int j=0; j<bundles[i].projections.size(); j++)
		{
			m = m + bundles[i].projections[j];
		}
	
		m = m/points.size();
		
		f32 s = 0.0;
		
		for(int j=0; j<bundles[i].projections.size(); j++)
			s += pow((bundles[i].projections[j][0][0] - m[0][0]),2) + 
					 pow((bundles[i].projections[j][1][0] - m[1][0]),2);
	
		s = sqrt(s/(2*bundles[i].projections.size()));

		Matrix T(3,3, {
			1/s, 0, -m[0][0]*1/s,
			0, 1/s, -m[1][0]*1/s,
			0, 0, 1
		});
	
		for(int j=0; j<bundles[i].projections.size(); j++)
		{
			bundles[i].projections[j] = T*bundles[i].projections[j];
		}
	}

	Matrix p11 = bundles[0].projections[0];
	Matrix p21 = bundles[0].projections[1];
	Matrix p31 = bundles[0].projections[2];
	Matrix p41 = bundles[0].projections[3];
	Matrix p51 = bundles[0].projections[4];
	Matrix p61 = bundles[0].projections[5];
	Matrix p71 = bundles[0].projections[6];
	Matrix p81 = bundles[0].projections[7];

	Matrix p12 = bundles[1].projections[0];
	Matrix p22 = bundles[1].projections[1];
	Matrix p32 = bundles[1].projections[2];
	Matrix p42 = bundles[1].projections[3];
	Matrix p52 = bundles[1].projections[4];
	Matrix p62 = bundles[1].projections[5];
	Matrix p72 = bundles[1].projections[6];
	Matrix p82 = bundles[1].projections[7];

	Matrix p1[8] = { p11, p21, p31, p41, p51, p61, p71, p81 };
	Matrix p2[8] = { p12, p22, p32, p42, p52, p62, p72, p82 };

	eightPointAlgorithm(p1, p2);
	
	return 0;
}
