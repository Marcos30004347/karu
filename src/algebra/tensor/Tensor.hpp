
typedef float f32;
typedef unsigned u32;

#include <cmath>
#include <assert.h>
#include <iostream>

#define min(x, y) x < y ? x : y
#define max(x, y) x < y ? y : x

class Tensor
{
	public:

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

	Tensor(u32 rank, const u32* dims);
	Tensor(std::initializer_list<u32> dims);
};

Tensor::Tensor(u32 rank, const u32* dims)
{
	this->rank = rank;
	this->dims = new u32[this->rank];
	this->stride = new u32[this->rank];
	
	for(int i=0; i<this->rank; i++)
	{
		this->dims[i] = dims[i];
	}

	this->page_h = 4;
	this->page_w = 4;

	this->srows = ceil(this->dims[this->rank - 2]/(f32)this->page_h) * this->page_h;
	this->scols = ceil(this->dims[this->rank - 1]/(f32)this->page_w) * this->page_w;


	if(rank >= 1)
	{
		this->stride[rank - 1] = 1;
	}

	if(rank >= 2)
	{
		this->stride[rank - 2] = this->scols;
	}

	if(rank >= 3)
	{
		this->stride[rank - 3] = this->stride[rank - 2] * this->srows;
	}

	for(int i = rank - 4; i >= 0; i--)
	{
		this->stride[i] = this->stride[i - 1] * this->dims[i];
	}

	this->size = this->stride[0] * this->dims[0];
	this->data = new f32[this->size];

	for(int i=0; i<this->size; i++)
	{
		this->data[i] = 0.0;
	}
}

Tensor::Tensor(std::initializer_list<u32> dims)
{
	this->rank = dims.size();
	this->dims = new u32[this->rank];
	this->stride = new u32[this->rank];
	
	for(int i=0; i<this->rank; i++)
	{
		this->dims[i] = dims.begin()[i];
	}

	this->page_h = 4;
	this->page_w = 4;

	this->srows = ceil(this->dims[this->rank - 2]/(f32)this->page_h) * this->page_h;
	this->scols = ceil(this->dims[this->rank - 1]/(f32)this->page_w) * this->page_w;


	if(rank >= 1)
	{
		this->stride[rank - 1] = 1;
	}

	if(rank >= 2)
	{
		this->stride[rank - 2] = this->scols;
	}

	if(rank >= 3)
	{
		this->stride[rank - 3] = this->stride[rank - 2] * this->srows;
	}

	for(int i = rank - 4; i >= 0; i--)
	{
		this->stride[i] = this->stride[i - 1] * this->dims[i];
	}

	this->size = this->stride[0] * this->dims[0];
	this->data = new f32[this->size];

	for(int i=0; i<this->size; i++)
	{
		this->data[i] = 0.0;
	}
}

f32 rec_get(f32* page, u32* stride, const u32* idx, u32 ph, u32 pw, u32 sr, u32 sc,  u32 rank, u32 j)
{
	// Block Matrix access
	if(j == rank - 2)
	{
		int by = idx[rank - 2] / ph;
		int bx = idx[rank - 1] / pw;
		int y = idx[rank - 2] - by * ph;
		int x = idx[rank - 1] - bx * pw;
		int bl = sc / pw;
		int bs = ph * pw;

		return page[(by*bl*bs + bx*bs) + y * pw + x];
	}

	// Block Array access
	if(j == rank - 1)
	{
		int bx = idx[rank - 1] / pw;
		int x = idx[rank - 1] - bx * pw;
		int bl = sc / pw;
		int bs = ph * pw;

		return page[bx*bs + x];
	}

	return rec_get(&page[idx[j] * stride[j]], stride, idx, ph, pw, sr, sc, rank, j + 1);
}

f32 get(Tensor t, u32* idx)
{
	return rec_get(t.data, t.stride, idx, t.page_h, t.page_w, t.srows, t.scols, t.rank, 0);
}

f32 get(Tensor t, std::initializer_list<u32> idx)
{
	return rec_get(t.data, t.stride, idx.begin(), t.page_h, t.page_w, t.srows, t.scols, t.rank, 0);
}

f32 rec_set(f32* page, f32 val, u32* stride, const u32* idx, u32 ph, u32 pw, u32 sr, u32 sc,  u32 rank, u32 j)
{
	// Block Matrix access
	if(j == rank - 2)
	{
		int by = idx[rank - 2] / ph;
		int bx = idx[rank - 1] / pw;
		int y = idx[rank - 2] - by * ph;
		int x = idx[rank - 1] - bx * pw;
		int bl = sc / pw;
		int bs = ph * pw;

		return page[(by*bl*bs + bx*bs) + y * pw + x] = val;
	}

	// Block Array access
	if(j == rank - 1)
	{
		int bx = idx[rank - 1] / pw;
		int x = idx[rank - 1] - bx * pw;
		int bl = sc / pw;
		int bs = ph * pw;

		return page[bx*bs + x] = val;
	}

	return rec_set(&page[idx[j] * stride[j]], val, stride, idx, ph, pw, sr, sc, rank, j + 1);
}

f32 set(Tensor t, u32* idx, f32 val)
{
	return rec_set(t.data, val, t.stride, idx, t.page_h, t.page_w, t.srows, t.scols, t.rank, 0);
}

f32 set(Tensor t, std::initializer_list<u32> idx, f32 val)
{
	return rec_set(t.data, val, t.stride, idx.begin(), t.page_h, t.page_w, t.srows, t.scols, t.rank, 0);
}


// Map a global idx to n tensor indices
void mapIdToTensorIdx(u32 idx, u32* dims, u32* t_idx, u32 n)
{
	int* cache = new int[n + 1];

	cache[n] = idx;

	// Retrieve indices of the b tensor
	for(int i = n - 1; i >= 0; i--)
	{
		cache[i] = cache[i + 1] / dims[i];
		t_idx[i] = cache[i + 1] % dims[i];
	}

	delete[] cache;
}

Tensor dot(Tensor a, Tensor b, u32 n, u32* a_idx, u32* b_idx)
{
	for(int i=0; i < n; i++)
	{
		assert(a.dims[a_idx[i]] == b.dims[b_idx[i]]);
	}

	// u[i] is strue if the i dimensions is being reduced in the A Tensor
	bool* u = new bool[a.rank];
	// initialize all values as false
	for(int i=0; i < a.rank; i++)
	{
		u[i] = false;
	}
	
	// p[i] is strue if the i dimensions is being reduced in the B Tensor
	bool* p = new bool[b.rank];
	// initialize all values as false
	for(int i=0; i < b.rank; i++)
	{
		p[i] = false;
	}


	// aidx[i] will store the current idx of of the A Tensor in the i'th dimension
	int* aidx = new int[a.rank];

	// bidx[i] will store the current idx of of the B Tensor in the i'th dimension
	int* bidx = new int[b.rank];

	// idxs[i] will store the current idx of the i'th dimension being reduced
	int* idxs = new int[n];
	
	// Helper array that is gong to hold the size
	// of the dimensions being reduced, should
	// be indexed by idx like dim[idx[i]], because
	// the dimensions being reduced have the same 
	// size in the A Tensor and in the B Tensor,
	// they can be picked by both by a_idx and a.dims
	// and be b_idx and b.dims
	u32* dim = a.dims;
	u32* idx = a_idx;


	/**
	/* Step 1: compute how many reductions will be needed
	**/
	for(int i=0; i<n; i++)
	{
		u[a_idx[i]] = true;
		p[b_idx[i]] = true;
	}


	// Number of accesses in the A tensor
	int a_ops = 1;
	for(int i=0; i < a.rank; i++)
	{
		a_ops *= 1*u[i] + !u[i]*a.dims[i];
	}

	// Number of accesses in the B tensor
	int b_ops = 1;
	for(int i=0; i < b.rank; i++)
	{
		b_ops *= 1*p[i] + !p[i]*b.dims[i];
	}

	// Number of writes in the C tensor
	int ops = a_ops * b_ops;

	/**
	/* Stage 2: Setup the resulting tensor
	**/
	unsigned* c_dims = new unsigned[a.rank - n + b.rank - n];
	
	// Setup the C tensor dimensions
	int q = 0;
	
	for(int i=0; i < a.rank; i++)
	{
		while(u[i])
		{
			q++;
			i++;
		}
	
		c_dims[i - q] = a.dims[i];
	}
	
	q = 0;
	
	for(int i=0; i < b.rank; i++)
	{
		while(p[i])
		{
			q++;
			i++;
		}

		c_dims[a.rank - n + i - q] = b.dims[i];
	}

	// Create the resulting Tensor 
	Tensor c(a.rank - n + b.rank - n, c_dims);

	/**
	/* Step 3: perform the reduction
	**/
	// The following loop can be parallelized
	for(int id = 0; id < ops; id++)
	{
		u32* cidx = new u32[c.rank];
		
		// Convert the loop id into the C Tensor indexes
		mapIdToTensorIdx(id, c.dims, cidx, c.rank);

		for(int i=0; i<n; i++)
		{
			idxs[i] = 0;
		}
	
		// TODO: REMOVE DEBUG PRINT
		{
			std::cout << "c";
			for(int i= 0; i<c.rank; i++)
			{
				std::cout << "[" << cidx[i] << "]";
			}
			std::cout << " = ";
		}

		for(;;)
		{

			int j = 0;
			// Get the A Tensor and B Tensor indexes
			for(int i=0; i < a.rank; i++)
			{
				if(u[i])
				{
					aidx[i] = idxs[j];
					bidx[i] = idxs[j];
					j = j + 1;
				}
				else
				{
					aidx[i] = cidx[i];
					bidx[i] = cidx[i + a.rank - n];
				}
			}

			// TODO: REMOVE DEBUG PRINT
			{
				std::cout << "a";
				for(int i=0; i<a.rank; i++)
				{
					std::cout << "[" << aidx[i] << "]";
				}
				std::cout << " * ";
				std::cout << "b";
				for(int i=0; i<b.rank; i++)
				{
					std::cout << "[" << bidx[i] << "]";
				}
				std::cout << "\n";
			}


			// Compute the number of dimensions
			// that where fully reduced
			int t = 0;	
			for(int i = 0; i < n; i ++)
			{
				if(idxs[i] == dim[idx[i]] - 1)
				{
					t += 1;
				}
			}

			// if all dimensions where fully reduced break
			if(t == n)
			{
				break;
			}


			// TODO: REMOVE DEBUG PRINT
			{
				std::cout << "		 +";
			}

			// Increase the index of the reduced dimensions
			idxs[n - 1] += 1;

			for(int i = n - 1; i >= 0; i --)
			{
				if(idxs[i] == dim[idx[i]])
				{
					idxs[i - 1] += 1;
					idxs[i] = 0;
				}
			}

		}
	
		// TODO: REMOVE DEBUG PRINT
		{
			std::cout << "\n";
		}
	
		delete []cidx;
	}

	delete[] u;
	delete[] p;
	delete[] aidx;
	delete[] bidx;
	delete[] idxs;

	return c;
}
