#include "algebra/matrix/MatrixEchelonForm.hpp"

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

void MatrixEchelonForm::toEchelonForm(MatrixData* M)
{
	u64 lead = 0;
	u64 rowCount = M->lines();
	u64 columnCount = M->columns();

	for(i64 r=0; r<rowCount; r++)
	{
		if(columnCount <= lead) break;
		
		i64 i = r;
		
		while(M->get(i, lead) == 0)
		{
			i = i + 1;
			if(rowCount == i)
			{
				i = r;
				lead = lead + 1;
				if(columnCount == lead)
					break;
			}
		}

		swapLines(M, i, r);
		
		u64 tmp = i;
		
		if(M->get(r, lead) != 0)
		{
			f32 m = M->get(r,lead);
			for(i32 j=0; j<columnCount; j++)
				M->set(r,j, M->get(r,j)/m);
		}

		for(i64 i=0; i<rowCount; i++)
		{
			if(i != r)
			{
				f32 m = M->get(i, lead);
				for(i32 j=0; j<columnCount; j++)
					M->set(i,j, M->get(i,j) - M->get(r,j)*m);
			}
		}
	
		lead = lead+1;
	}
}

}
