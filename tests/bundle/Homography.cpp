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

	std::vector<Matrix> projections[2];

	for(u64 i=0; i<2; i++)
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

		for(int j=0; j<points.size(); j++)
		{
			Matrix projection = camera.projection(points[j]);
			projections[i].push_back(projection);
			point_idx.push_back(j);
		}
	}

	Matrix p1[8] = {
		projections[0][0],
		projections[0][1],
		projections[0][2],
		projections[0][3],
		projections[0][4],
		projections[0][5],
		projections[0][6],
		projections[0][7]
	};

	Matrix p2[8] = {
		projections[1][0],
		projections[1][1],
		projections[1][2],
		projections[1][3],
		projections[1][4],
		projections[1][5],
		projections[1][6],
		projections[1][7]
	};

	Matrix F = eightPointAlgorithm(p1, p2);

	f32 f1, f2;

	f32 a1 = 1000/1000;
	f32 a2 = 1000/1000;

	estimateCameraFocalLengths(F, Matrix(3,1, { 500, 500, 1 }), Matrix(3,1, { 500, 500, 1 }), &f1, &f2);

	Matrix K1(3,3,
	{
		f1*a1, 0, 500,
		0,    f1, 500,
		0,     0,   1
	});

	Matrix K2(3,3,
	{
		f2*a2, 0, 500,
		0,    f2, 500,
		0,     0,   1
	});

	Matrix E = getEssentialMatrix(F, K1, K2);

	Matrix R1, R2, t1, t2;

	estimateRotationAndTranslation(E, R1, R2, t1, t2);

	Matrix R[2] = {R1, R2};
	Matrix T[2] = {t1, t2};
	Matrix I[2] = {K1, K2};

	Matrix points1[8] = {
		projections[0][0],
		projections[0][1],
		projections[0][2],
		projections[0][3],
		projections[0][4],
		projections[0][5],
		projections[0][6],
		projections[0][7]
	};

	Matrix points2[8] = {
		projections[1][0],
		projections[1][1],
		projections[1][2],
		projections[1][3],
		projections[1][4],
		projections[1][5],
		projections[1][6],
		projections[1][7]
	};

	chooseRealizableSolution(I, R, T, points1, points2, 8);

	return 0;
}
