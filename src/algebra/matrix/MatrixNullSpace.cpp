#include <vector>

#include "algebra/matrix/MatrixNullSpace.hpp"
#include "algebra/matrix/MatrixEchelonForm.hpp"
#include "algebra/matrix/Matrix.hpp"

namespace karu::algebra  {

void MatrixNullSpace::nullSpace(MatrixData* M, MatrixData& ns) {
	
	MatrixEchelonForm::toEchelonForm(M);
	std::cout << M->lines() << " " << M->columns() << "\n";
	Matrix _m(*M);
	printMatrix(_m);
	std::cout << _m.rows() << " " << _m.columns() << "\n";

	u64 columns = M->columns();
	u64 lead = 0;

	while(M->get(lead, lead) == 1 && lead < M->lines())
		lead++;

	u64 rank = columns - lead;

	if(rank == 0)
	{
		std::vector<f32> zeros(columns, 0);
		ns = MatrixData(1,columns,1,1, zeros.data());
		return;
	}

	std::vector<f32> basis(rank*columns, 0);
	// std::cout << "\n";
	// std::cout << "columns " << columns << "\n";
	// std::cout << "rank " << rank << "\n";

	for(i64 i=1; i<=rank; i++)
		basis[i*columns - rank + (i-1)] = 1;
	
	for(i64 i=0; i<rank; i++)
		for(i64 j=0; j<M->lines(); j++)
			basis[i*columns + j] -= M->get(j, lead+i);
		
	ns = MatrixData(rank, columns, 1, 1, basis.data());
}

}
