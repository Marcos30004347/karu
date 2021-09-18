
typedef float f32;
typedef unsigned u32;

#include <cmath>

#define min(x, y) x < y ? x : y
#define max(x, y) x < y ? y : x

struct Tensor
{
	// Tensor data
	f32* data;
	// Tensor rank
	u32  rank;
	// Stride is an array that contains
	// how many elements should be skiped
	// to get the starting point of the 
	// elements one unit in that dimentsion
	//
	// For instance, given a tensor of rank 3 T_ikj,
	// stride[1] will tell us how many elements
	// should be skiped from T_ikj to get to T_i(j+1)k
	u32* stride;
	
	// Tensor dimension sizes
	u32* dims;

	// Ammount of elements stored inside the Tensor
	u32  size;

	// If the tensor have two or more dimensions,
	// the last two can be seen as a matrix. 

	// The Width and Height of blocks of the last two
	// dimensions Matrices
	u32  page_w;
	u32  page_h;
	
	//srows and scols store how many rows anc columns are
	// actualy stored on those matrices.
  // Becuase the Tensor is stored as matrices blocks, srows 
	// and scols can differ from dims[rank - 2] and dims[rank - 1] 
	u32  srows;
	u32  scols;

};

Tensor createTensor(u32 rank, u32* dims)
{
	Tensor t;

	t.rank = rank;
	t.dims = new u32[t.rank];
	t.stride = new u32[t.rank];

	t.page_h = 16;
	t.page_w = 16;

	t.srows = ceil(t.dims[t.rank - 2]/(f32)t.page_h) * t.page_h;
	t.scols = ceil(t.dims[t.rank - 1]/(f32)t.page_w) * t.page_w;

	t.stride[0] = dims[0];

	t.size = dims[0];
 
	for(u32 i=1; i < rank - 2; i++)
	{
		t.size *= dims[i];
	}

	t.size *= t.srows;
	t.size *= t.scols;


	if(rank >= 1)
	{
		t.stride[rank - 1] = t.scols;
	}

	if(rank >= 2)
	{
		t.stride[rank - 2] = t.stride[rank - 1] * t.srows;
	}

	for(u32 i=rank - 3; i >=0 ; i--)
	{
		t.stride[i] = t.stride[i - 1] * t.dims[i];
	}

	for(u32 i=1; i<rank - 2; i++)
	{
		t.stride[i] += t.stride[i - 1] + t.dims[i];
	}

	for(u32 i=0; i<rank; i++)
	{
		t.dims[i] = t.dims[i];
	}

	t.data = new f32[t.size];

	for(u32 i=0; i<t.size; i++)
	{
		t.data[i] = 0.0;
	}


	return t;
}


f32 get(Tensor t, u32* idx)
{
	u32 rank = t.rank;

	f32* page = t.data;

	if(rank == 1)
	{
		return page[idx[0]];
	}

	if(rank > 2)
	{
		u32 stride;

		// discover how many elements to skip
		// to get to the last two dimensions
		// that can be seen as a matrix
		for(u32 d = 0; d < rank - 2; d++)
		{
			for(u32 i = 0; i < idx[d]; i++)
			{
				stride += t.stride[d];
			}
		}

		page = &t.data[stride];
	}

	u32 by = idx[rank - 2]/t.page_h;
	u32 bx = idx[rank - 1]/t.page_w;

	u32 y = idx[rank - 2] - by * t.page_h;
	u32 x = idx[rank - 1] - bx * t.page_w;

	u32 bl = t.scols / t.page_w;
	u32 bs = t.page_h * t.page_w;

	u32 stride = by*bl*bs + bx*bs;

	return page[stride + y * t.page_w + x];
}

f32 set(Tensor t, u32 rank, u32* idx, f32 val)
{
	f32* page = t.data;

	if(rank == 1)
	{
		return page[idx[0]];
	}

	if(rank > 2)
	{
		u32 stride;

		// discover how many elements to skip
		// to get to the last two dimensions
		// that can be seen as a matrix
		for(u32 d = 0; d < rank - 2; d++)
		{
			for(u32 i = 0; i < idx[d]; i++)
			{
				stride += t.stride[d];
			}
		}

		page = &t.data[stride];
	}

	u32 by = idx[rank - 2]/t.page_h;
	u32 bx = idx[rank - 1]/t.page_w;

	u32 y = idx[rank - 2] - by * t.page_h;
	u32 x = idx[rank - 1] - bx * t.page_w;

	u32 bl = t.scols / t.page_w;
	u32 bs = t.page_h * t.page_w;

	u32 stride = by*bl*bs + bx*bs;

	return page[stride + y * t.page_w + x] = val;
}

Tensor dot(Tensor a, Tensor b, u32 n, u32* a_idx, u32* b_idx)
{
	// c_[i][j][k][e][f][q] = a_[i][j][k][t] * b_[e][f][t][q]
	// c_[i][j][k][e][f][q] = sum t a_[i][j][k][t]*b_[e][f][t][q]

	bool* a_red_vars = new bool[a.rank];
	bool* b_red_vars = new bool[b.rank];

	for(int i=0; i<a.rank; i++)
	{
		a_red_vars[i] = false;
	}

	for(int i=0; i<b.rank; i++)
	{
		b_red_vars[i] = false;
	}

	// Step 1: get how many reduction will be needed
	for(int i=0; i<n; i++)
	{
		a_red_vars[a_idx[i]] = true;
		b_red_vars[b_idx[i]] = true;
	}


	int a_ops = 1;
	for(int i=0; i < a.rank; i++)
	{
		if(!a_red_vars[i])
		{
			a_ops *= a.dims[i];
		}
	}

	int b_ops = 1;
	for(int i=0; i < b.rank; i++)
	{
		if(!b_red_vars[i])
		{
			b_ops *= b.dims[i];
		}
	}

	int ops = a_ops + b_ops;


	// Stage 2: Setup the resulting tensor
	unsigned* c_dims = new unsigned[a.rank - n + b.rank - n];
	
	int q;

	q = 0;
	for(int i=0; i < a.rank; i++)
	{
		while(a_red_vars[i])
		{
			q++;
		}

		c_dims[i - q] = a.dims[i];
	}

	q = 0;
	for(int i=0; i < b.rank; i++)
	{
		while(b_red_vars[i])
		{
			q++;
		}

		c_dims[a.rank - n + i - q] = b.dims[i];
	}

	Tensor c = createTensor(a.rank - n + b.rank - n, c_dims);


	// Step 3: perform the reduction
	// This stage can be parallelized
	int* tmp_a = new int[a.rank + 1];
	int* tmp_b = new int[b.rank + 1];

	for(int id = 0; id < ops; id++)
	{
		int* a_ijk = new int[a.rank];
		int* b_ijk = new int[b.rank];
	
		int g_id = id;
		int tops = ops;
		
		tmp_b[b.rank] = g_id;

		// Retrieve indices of the b tensor
		for(int i = b.rank - 1; i >= 0; i--)
		{
			int old_i = i;

			// Ignore variables that are being reduced
			while(b_red_vars[i])
			{
				b_ijk[i] = -1;
				tmp_b[i] = tmp_b[i + 1];
				i = i - 1;
			}
		
			tmp_b[i] = floor(tmp_b[i + 1] / b.dims[i]);
			b_ijk[i] = tmp_b[i + 1] % b.dims[i];
		}

		// Retrieve indices of the a tensor
		for(int i = a.rank - 1; i >= 0; i--)
		{
			int old_i = i;

			// Ignore variables that are being reduced
			while(a_red_vars[i])
			{
				a_ijk[i] = -1;
				tmp_a[i] = tmp_a[i + 1];
				i = i - 1;
			}

			tmp_a[i] = floor(tmp_a[i + 1] / a.dims[i]);
			a_ijk[i] = tmp_a[i + 1] % a.dims[i];
		}

		// a_ijk contain the index in every dimension
		// that are not being reduced of a.
		// b_ijk contain the index in every dimension
		// that are not being reduced of b.
		// Then c can be indexed by a[i] with i not existent in
		// a_idx and by b[i] with i not existent in b_idx
		// c[a_ijk[0]][a_ijk[1]]...[a_ijk[a.rank - n][b_ijk[0]][b_ijk[1]]...[b_ijk[b.rank - n]

		// TODO: Reduce tensor
	}

	// // Iterate over tensor a variables
	// for(int i=0; i < a.dims[0]; i++)
	// {
	// 	for(int j=0; j < a.dims[1]; j++)
	// 	{
	// 		for(int k=0; k < a.dims[2]; k++)
	// 		{

	// 			// Iterate over tensor b variables
	// 			for(int e=0; e < b.dims[0]; e++)
	// 			{
	// 				for(int f=0; f < b.dims[1]; f++)
	// 				{
	// 					for(int q=0; q < b.dims[3]; q++)
	// 					{

	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	return c;
}
