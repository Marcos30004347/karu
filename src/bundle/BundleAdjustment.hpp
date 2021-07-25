/**
** Bundle Adjustment Basics:

	problem:
		min f 
		where f = sum of i,m {(u[i][j] - proj(C[j], X[i]))**2}

		where u[i][j] is the observed point X[i] by the
		camera C[j].

	def:
		r = u[i][j] - proj(C[j], X[i])
	
		given x = [c, p] with 
		c being a flat vector with all the camera parameters and 
		p being a flat vector with all the point parameters. 

		lets define û[i][j] = proj(C[j], X[i])
		then define:
		 	A[i][j] = diff(û[i][j], c[j]), and diff(û[i][j], c[k]) = 0 for every j different from k 
 		 	B[i][j] = diff(û[i][j], p[i]), and diff(û[i][j], p[k]) = 0 for every i different from k 

	U[j] = sum over i for each point i that is seen by camera j of 	A[i][j].T * A[i][j]
	V[i] = sum over j for each camera j that sees the point i of 		B[i][j].T * B[i][j]
	W[i][j] = A[i][j].T * B[i][j]

	r[i][j] = u[i][j] - û[i][j]
	e[c[j]] = sum over all i seen in camera j of A[i][j]*r[i][j]
	e[p[i]] = sum over all j that sees point i of B[i][j]*r[i][j]

	Then we can get delta x by
	[[U W], [W.T, V]] * delta(x) = flat([e[c], e[p]])

	It is possibly to get 
	(U - W*inverse(V)*W.T) = S
	with being c is 1: Compute the derivation matrix:
Aij =
∂uˆij
∂cj
, Bij =
∂uˆij
∂Xi
rij = uij − π(cj , Xi) = uij − uˆij


/**
Algorithm: Compute Sx with W, V and x
without forming S on memory:
	x1 = W.T*x
	x2 = inverse(V) * x1
	x3 = W*x2
	x4 = B*x
	Sx = x4 - x3
**/

typedef int Vector;
typedef int Matrix;

float dot(Vector a, Vector b){}
float norm(Vector a){}

// solve A*x = b
Matrix ConjugateGradient(Matrix A, Matrix x0, Matrix b, float e)
{
	Matrix x_k, x_k_min_1, x_k_add_1;
	Matrix p_k, p_k_min_1, p_k_add_1;
	Matrix r_k, r_k_min_1, r_k_add_1;

	float alpha_k, alpha_k_add_1;
	float beta_k;

	x_k_min_1 = x0;
	r_k_min_1 = b - A*x_k_min_1;
	p_k = r_k_min_1;

	Matrix w = A*p_k;

 	alpha_k = dot(r_k_min_1, r_k_min_1)/(dot(p_k, w));

	x_k = x_k_min_1 + alpha_k*p_k;
	r_k = r_k_min_1 - alpha_k*w;

	int k = 1;

	while (norm(r_k) > e)
	{
		beta_k = dot(r_k, r_k)/(dot(r_k_min_1, r_k_min_1));
		p_k_add_1 = r_k + beta_k*p_k;
		
		w = A*p_k_add_1;
		
		alpha_k_add_1 = dot(r_k, r_k)/(dot(p_k_add_1, w));

		x_k_add_1 = x_k + alpha_k_add_1*p_k_add_1;
		r_k_add_1 = r_k - alpha_k_add_1*w;
		
		x_k_min_1 = x_k;
		x_k = x_k_add_1;
	
		r_k_min_1 = r_k;
		r_k = r_k_add_1;
	
		p_k_min_1 = p_k;
		p_k = p_k_add_1;
	
		alpha_k = alpha_k_add_1;
	}

	return x_k;
}


// solve A*x0 = b
void PreconditionedConjugateGradient(
	Matrix M,
	Matrix A,
	Vector x0,
	Vector b
)
{
	Vector r0 = b - A*x0;

}
