#include "MatrixLU.hpp"

#include <math.h> 

namespace karu::algebra
{

void MatrixLU::LUdecompose(MatrixData* L, MatrixData* U, const MatrixData* const A)
{
	i32 i = 0, j = 0, k = 0;

	i32 n = A->lines();

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (j < i)
			{
				L->set(j, i, 0);
			}
			else
			{
				L->set(j, i, A->get(j, i));
				for (k = 0; k < i; k++) {
					L->set(j,i, L->get(j,i) - L->get(j,k)*U->get(k,i));
				}
			}
		}
		for (j = 0; j < n; j++)
		{
				if (j < i)
				{
					U->set(i, j, 0);
				}
				else if (j == i)
				{
					U->set(i,j,1); 
				}
				else {
					U->set(i,j, A->get(i,j)/L->get(i,i));

					for (k = 0; k < i; k++) {
							U->set(i,j, U->get(i,j) - ((L->get(i,k)*U->get(k,j))/L->get(i,i)));
					}
				}
		}
	}
}

void pivot(MatrixData* A, i32 i, i32 k)
{
	for(i32 j=0; j<A->columns(); j++)
	{
		f32 tmp = A->get(i, j);

		A->set(i,j, A->get(k, j));
		A->set(k,j, tmp);
	}
}

int MatrixLU::LUPdecompose(MatrixData* A, MatrixData* P)
{
	i32 i, j, k, imax;
	f32 maxA, absA;
	i32 N = A->lines();
	
	for(i=0; i<=N; i++)
	{
		P->set(i,0, i);
	}
	
	for(i=0; i<N; i++)
	{
		maxA = 0.f;
		imax = i;

		for(k=i;k<N;k++)
		{
			if((absA = fabs(A->get(k,i))) > maxA)
			{
				maxA = absA;
				imax = k;
			}
		}

		if(maxA < 0.00001) return 0;

		if(imax != i)
		{
			j = P->get(i,0);
			P->set(i, 0, P->get(imax, 0));
			P->set(imax, 0, j);
			pivot(A, i, imax);
			P->set(N, 0, P->get(N,0)+1);
		}

		for(j=i+1; j<N; j++)
		{
			A->set(j,i, A->get(j,i)/A->get(i,i));
			for(k=i+1; k<N; k++)
			{
				A->set(j,k, A->get(j,k) - A->get(j,i)* A->get(i,k));
			}
		}
	}
	return 1;
}

int MatrixLU::LUPSolve(const MatrixData* const A, const MatrixData* const P,  const MatrixData* const b,  MatrixData* x)
{
	i32 N = A->lines();

	for(i32 i=0; i<N; i++)
	{
		x->set(i,0, b->get(P->get(i,0), 0));
		for(i32 k=0; k<i; k++)
		{
			x->set(i,0, x->get(i,0) - A->get(i, k) * x->get(k,0));
		}
	}

	for(i32 i=N-1; i>=0; i--)
	{
		for(i32 k=i+1; k<N; k++)
		{
			x->set(i, 0, x->get(i,0) - A->get(i,k)*x->get(k,0));
		}
		x->set(i, 0, x->get(i,0)/A->get(i,i));
	}

	return 1;
}

i32 MatrixLU::LUPInvet(const MatrixData* const A, const MatrixData* const P,  MatrixData* A_Inv)
{
	i32 N = A->lines();

	for(i32 j=0; j<N; j++)
	{
		for(i32 i=0; i<N; i++)
		{
			A_Inv->set(i, j, static_cast<i32>(P->get(i,0)) == j ? 1.f : 0.f);

			for(i32 k=0; k<i; k++)
			{
				A_Inv->set(i, j, A_Inv->get(i,j) - A->get(i,k)*A_Inv->get(k,j));
			}
		}

		for(i32 i=N-1; i>=0; i--)
		{
			for(i32 k=i+1; k<N; k++)
			{
				A_Inv->set(i, j, A_Inv->get(i,j) - A->get(i,k)*A_Inv->get(k,j));
			}
			A_Inv->set(i, j, A_Inv->get(i,j) / A->get(i,i));

		}
	}

	return 1;
}

f32 MatrixLU::LUPDeterminant(const MatrixData* const A, const MatrixData* const P)
{
	f32 det = A->get(0,0);
	i32 N = A->lines();

	for(i32 i=1; i<N; i++)
	{
		det *= A->get(i, i);
	}

	return (static_cast<i32>(P->get(N,0)) - N)%2 == 0 ? det : -det;
}

}
