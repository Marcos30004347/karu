#include "algebra/matrix/MatrixEchelonForm.hpp"
#include <iomanip>
#include <cmath>
namespace karu::algebra {

void swapLines(MatrixData* A, i32 i, i32 k)
{
	for(i32 j=0; j<A->columns(); j++)
	{
		f32 tmp = A->get(i, j);

		A->set(i,j, A->get(k, j));
		A->set(k,j, tmp);
	}
}

u64 findPivo(MatrixData* M, u64 l, u64 c)
{
	u64 best = l;
	f32 val = fabs(M->get(l,c));

	for(u64 r=l+1; r<M->lines(); r++)
	{
		f32 curr = fabs(M->get(r, c));
		if(curr > val)
		{
			best = r;
			val = curr;
		}
	}

	return best;
}

void MatrixEchelonForm::toEchelonForm(MatrixData* M)
{
	u64 rows = M->lines();
	u64 cols = M->columns();

	u64 r = 0;

	for(u64 c=0; c<cols; c++)
	{

		u64 pivo = findPivo(M, r, c);
		if(pivo != r) {
			// std::cout << "swap " << r << " " << pivo << "\n";
			swapLines(M, r, pivo);
			// std::cout << "now\n";
			// for(i64 i=0; i<rows; i++)
			// {
			// 	for(i64 j=0; j<cols; j++)
			// 	{
			// 		std::cout << std::defaultfloat << M->get(i,j) << " ";
			// 	}
			// 	std::cout << "\n";
			// }
			// std::cout << "\n";
			// for(i64 i=0; i<rows; i++)
			// {
			// 	for(i64 j=0; j<cols; j++)
			// 	{
			// 		std::cout << std::defaultfloat << M->get(i,j) << " ";
			// 	}
			// 	std::cout << "\n";
			// }
			// std::cout << "\n";
		}

		f32 b = M->get(r, c);
		
		for(u64 j=c; j<cols; j++)
		{
			f32 a = M->get(r, j);
			M->set(r, j, a/b);
		}

		if(r > 0)
		{
			for(u64 l=0; l<r; l++)
			{
				f32 k = M->get(l, c);
				for(u64 j=c; j<cols; j++)
				{
					f32 a = M->get(l, j);
					f32 b = M->get(r, j);

					M->set(l, j, a - b*k);
				}
			}
			// std::cout << "above\n";
			// for(i64 i=0; i<rows; i++)
			// {
			// 	for(i64 j=0; j<cols; j++)
			// 	{
			// 		std::cout << std::defaultfloat << M->get(i,j) << " ";
			// 	}
			// 	std::cout << "\n";
			// }
			// std::cout << "\n";
		}


		if(r < rows-1)
		{

			for(u64 l=r+1; l<rows; l++)
			{
				f32 k = M->get(l, c);
				for(u64 j=c; j<cols; j++)
				{
					f32 a = M->get(l, j);
					f32 b = M->get(r, j);

					M->set(l, j, a - b*k);
				}
			}
			// std::cout << "bellow\n";
			// for(i64 i=0; i<rows; i++)
			// {
			// 	for(i64 j=0; j<cols; j++)
			// 	{
			// 		std::cout << std::defaultfloat << M->get(i,j) << " ";
			// 	}
		// 	// 	std::cout << "\n";
		// 	// }
		// 	// std::cout << "\n";
		}

		// for(i64 i=0; i<rows; i++)
		// {
		// 	for(i64 j=0; j<cols; j++)
		// 	{
		// 		std::cout << std::defaultfloat << M->get(i,j) << " ";
		// 	}
		// 	std::cout << "\n";
		// }
		// std::cout << "\n";
		// std::cout << "\n";
		r+=1;

		if(r == rows)
			break;
	}
	// u64 lead = 0;
	// u64 rowCount = M->lines();
	// u64 columnCount = M->columns();
	// std::cout << "\n";
	// for(i64 i=0; i<rowCount; i++)
	// {
	// 	for(i64 j=0; j<columnCount; j++)
	// 	{
	// 		std::cout << std::fixed << std::setprecision(30) << M->get(i,j) << " ";
	// 	}
	// 	std::cout << "\n";
	// }
	// std::cout << "\n";
	// for(i64 r=0; r<rowCount; r++)
	// {
	// 	if(columnCount <= lead) break;
		
	// 	i64 i = r;
		
	// 	while(M->get(i, lead) == 0)
	// 	{
	// 		i = i + 1;
	// 		if(rowCount == i)
	// 		{
	// 			i = r;
	// 			lead = lead + 1;
	// 			if(columnCount == lead)
	// 				break;
	// 		}
	// 	}
	
	// 	if(i != r) swapLines(M, i, r);
		
	// 	f32 tmp = M->get(r, lead);
	
	// 	if(tmp != 0.0)
	// 	{
	// 		for(i32 j=0; j<columnCount; j++)
	// 			M->set(r,j, M->get(r,j)/tmp);
	// 	}

	// 	for(i64 j=0; j<columnCount; j++)
	// 	{
	// 		std::cout << M->get(r,j) << " ";
	// 	}
	// 	std::cout << "\n";
	// 	std::cout << "\n";

	// 	for(i64 i=0; i<rowCount; i++)
	// 	{
	// 		if(i != r)
	// 		{
	// 			f32 k = M->get(i, lead);
	// 			for(i32 j=0; j<columnCount; j++)
	// 			{
	// 				// std::cout << M->get(i,j) << " - " << M->get(r,j) << " * " << k << " = " << M->get(i,j) - M->get(r,j)*k <<  "\n";

	// 				M->set(i,j, M->get(i,j) - M->get(r,j)*k);
	// 			}
	// 		}
	// 	}

	// 	for(i64 i=0; i<rowCount; i++)
	// 	{
	// 		for(i64 j=0; j<columnCount; j++)
	// 		{
	// 			std::cout << std::defaultfloat << M->get(i,j) << " ";
	// 		}
	// 		std::cout << "\n";
	// 	}
	// 	std::cout << "\n";

	// 	lead = lead+1;
	// }
}

}
