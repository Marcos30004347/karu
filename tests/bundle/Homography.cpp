#include <assert.h>

#include "bundle/Homography.hpp"
#include "bundle/BundleAdjustment.hpp"
#include "renderer/Renderer.hpp"
using namespace karu;
using namespace karu::bundle;
using namespace karu::algebra;


int main()
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

	std::vector<Bundle> bundles;

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


	for(u64 i=0; i<4; i++)
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
		std::vector<Matrix> projections;
	
		// std::vector<Matrix> pixels;

		for(int j=0; j<points.size(); j++)
		{
			Matrix projection = camera.projection(points[j]);
			// pixels.push_back(Matrix(2,1, {projection[0][0], projection[1][0]}));
			point_idx.push_back(j);
			projections.push_back(projection);
		}
		// Renderer r(1000,1000);
		// r.draw2dPoints(camera, pixels);

		bundles.push_back({camera, projections, point_idx});
	}

	Matrix K(3,3, 
	{
		bundles[0].camera.fx, 0, bundles[0].camera.cx,
		0, bundles[0].camera.fy, bundles[0].camera.cy,
		0, 										0, 										1
	});


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
	Matrix p1_[8] = { p11, p21, p31, p41, p51, p61, p71, p81 };
	Matrix p2_[8] = { p12, p22, p32, p42, p52, p62, p72, p82 };

	Matrix F = eightPointAlgorithm(p1, p2);

	std::cout << "Fundamental:\n";
	printMatrix(F);

	std::cout << "errors Fundamental:\n";
	for(i64 i=0;i<8; i++)
		printMatrix(transpose(p2_[i])*F*p1_[i]);

	printf("=============\n");
	estimateCameraFocalLengths(F);
	printf("=============\n");

	Matrix E = getEssentialMatrix(F, K, K);

	printMatrix(E);

	Matrix R1, R2, t1, t2;

	estimateRotationAndTranslation(E, R1, R2, t1, t2);

	f32 x, y, z;

	rotationToEulerAngles(R1, &x, &y, &z);
	rotationToEulerAngles(R2, &x, &y, &z);
	std::cout << "t1\n";
	printMatrix(t1);
	std::cout << "t2\n";
	printMatrix(t2);
	std::cout << "\n";

	Matrix R[2] = {R1, R2};
	Matrix T[2] = {t1, t2};
	Matrix I[2] = {K, K};
	
	chooseRealizableSolution(I, R, T, p1_, p2_, 8);
	
	return 0;
}
