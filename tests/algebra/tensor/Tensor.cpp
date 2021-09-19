#include "algebra/tensor/Tensor.hpp"

int main()
{
	unsigned adims[3] = {3, 3, 3};
	unsigned bdims[3] = {3, 3, 3};

	Tensor a({3, 3, 3});
	Tensor b({3, 3, 3});

	set(a, { 0, 1, 0 }, 1.0);
	set(a, { 0, 1, 1 }, 1.0);
	set(a, { 1, 1, 1 }, 1.0);

	std::cout << get(a, {0, 1, 0}) << "\n";
	std::cout << get(a, {0, 1, 1}) << "\n";
	std::cout << get(a, {1, 1, 1}) << "\n";

	for(int i=0; i<a.size; i++)
	{
		std::cout << a.data[i] << " ";
	}
	std::cout << "\n";


	unsigned aidx[2] = {1, 2};
	unsigned bidx[2] = {1, 2};

	dot(a, b, 2, aidx, bidx);
	
	return 0;
}
