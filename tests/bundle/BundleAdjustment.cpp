#include <assert.h>

// #include "helpers.hpp"

#include "camera/Camera.hpp"
#include "renderer/Renderer.hpp"
#include "bundle/BundleAdjustment.hpp"
#include "algebra/linear/Rotation.hpp"
#include "algebra/compute/Compute.hpp"

using namespace karu;
using namespace karu::algebra;
using namespace karu::algebra::compute;
using namespace karu::bundle;


int main()
{
	// algebra::compute::Context::initContext();

	std::vector<Matrix> points;

	points.push_back(Matrix(3,1,{ 1.0,1.0,0.0 }));
	points.push_back(Matrix(3,1,{ 1.0,-1.0,0.0 }));
	points.push_back(Matrix(3,1,{ -1.0,-1.0,0.0 }));
	points.push_back(Matrix(3,1,{ -1.0,1.0,0.0 }));
	points.push_back(Matrix(3,1,{ 0.0, 0.0, sqrt(2.0) }));

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

		for(int j=0; j<5; j++)
		{
			Matrix projection = camera.projection(points[j]);

			point_idx.push_back(j);
			projections.push_back(projection);
		}

		Matrix position = positions[i] + position_noises[i];
		Matrix rotation = rotations[i] + rotation_noises[i];
		
		Camera cam(
			3000,
			3000,
			500,
			500, 
			position,
			rotation
		);
		
		bundles.push_back({cam, projections, point_idx});
	}

	SpMatrix U, V, V_inv, W, W_T;

	std::vector<std::vector<u64>> point_idx_to_camera(points.size());

	for(i64 j=0; j<bundles.size(); j++)
		for(u64 i : bundles[j].point_idx)
			point_idx_to_camera[i].push_back(j);


	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i=0; i<bundles[j].projections.size(); i++)
		{
				u64 k = bundles[j].point_idx[i];
				Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(points[k]);
				std::cout << norm(r) << "\n";
		}
	}

	f32 lambda = 0.0;

	do
	{
		solveNormalEquations(bundles, points, point_idx_to_camera, lambda);
		std::cout << "converging: " << lambda << std::endl;
	} while (lambda > 0.1);
	
	
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "\n";

	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i=0; i<bundles[j].projections.size(); i++)
		{
				u64 k = bundles[j].point_idx[i];
				Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(points[k]);
				// printMatrix(r);
				std::cout << "Reprojection Error: " << norm(r) << "\n";
		}
	}
	
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "\n";

	for(i64 j=0; j<bundles.size(); j++)
	{
		std::cout << "fx: " << bundles[j].camera.fx << "\n";
		std::cout << "fy: " << bundles[j].camera.fy << "\n";
		std::cout << "cx: " << bundles[j].camera.cx << "\n";
		std::cout << "cy: " << bundles[j].camera.cx << "\n";
		std::cout << "Cx: " << bundles[j].camera.P[0][0] << "\n";
		std::cout << "Cy: " << bundles[j].camera.P[1][0] << "\n";
		std::cout << "Cz: " << bundles[j].camera.P[2][0] << "\n";
		std::cout << "r1: " << bundles[j].camera.R[0][0] << "\n";
		std::cout << "r2: " << bundles[j].camera.R[1][0] << "\n";
		std::cout << "r3: " << bundles[j].camera.R[2][0] << "\n";
		std::cout << "\n";
	}

	// hessian(bundles, points, point_idx_to_camera, 0.0, U, V, W, W_T);
	// buildHessianVInverse(V, V_inv);
	
	// std::vector<u64> row_ptr = U.m_data.rowPtr();
	// for(int i=0; i<row_ptr.size(); i++)
	// {
	// 	std::cout << row_ptr[i] << " ";
	// }
	// std::cout << "\n";
	// std::vector<u64> col_idx = U.m_data.columnsIdx();
	// for(int i=0; i<col_idx.size(); i++)
	// {
	// 	std::cout << col_idx[i] << " ";
	// }
	// std::cout << "\n";
	// printMatrix(U);
	// std::cout << "\n";

	// row_ptr = W_T.m_data.rowPtr();
	// for(int i=0; i<row_ptr.size(); i++)
	// {
	// 	std::cout << row_ptr[i] << " ";
	// }
	// std::cout << "\n";
	// col_idx = W_T.m_data.columnsIdx();
	// for(int i=0; i<col_idx.size(); i++)
	// {
	// 	std::cout << col_idx[i] << " ";
	// }

	// std::cout << "\n";
	// printMatrix(W_T);
	// std::cout << "\n";

	// row_ptr = W.m_data.rowPtr();
	// for(int i=0; i<row_ptr.size(); i++)
	// {
	// 	std::cout << row_ptr[i] << " ";
	// }
	// std::cout << "\n";
	// col_idx = W.m_data.columnsIdx();
	// for(int i=0; i<col_idx.size(); i++)
	// {
	// 	std::cout << col_idx[i] << " ";
	// }

	// std::cout << "\n";
	// printMatrix(W);
	// std::cout << "\n";

	// row_ptr = V.m_data.rowPtr();
	// for(int i=0; i<row_ptr.size(); i++)
	// {
	// 	std::cout << row_ptr[i] << " ";
	// }
	// std::cout << "\n";
	// col_idx = V.m_data.columnsIdx();
	// for(int i=0; i<col_idx.size(); i++)
	// {
	// 	std::cout << col_idx[i] << " ";
	// }

	// std::cout << "\n";
	// printMatrix(V);
	// std::cout << "\n";

	// row_ptr = V_inv.m_data.rowPtr();
	// for(int i=0; i<row_ptr.size(); i++)
	// {
	// 	std::cout << row_ptr[i] << " ";
	// }
	// std::cout << "\n";
	// col_idx = V_inv.m_data.columnsIdx();
	// for(int i=0; i<col_idx.size(); i++)
	// {
	// 	std::cout << col_idx[i] << " ";
	// }
	// std::cout << "\n";
	// printMatrix(V_inv);
	// std::cout << "\n";

	// algebra::compute::Context::stopContext();
	return 0;
}
