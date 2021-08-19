#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sys/time.h>

#include "algebra/matrix/Matrix.hpp"

#define epsilon 1.e-8

using namespace std;
using namespace karu::algebra;

// #define SIGN(a,b) ((b) > 0.0 ? fabs(a) : - fabs(a))

// static double maxarg1,maxarg2;
// #define FMAX(a,b) (maxarg1 = (a),maxarg2 = (b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

// static int iminarg1,iminarg2;
// #define IMIN(a,b) (iminarg1 = (a),iminarg2 = (b),(iminarg1 < (iminarg2) ? (iminarg1) : iminarg2))

// static double sqrarg;
// #define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

// int svdcmp(double **a, int nRows, int nCols, double *w, double **v);

// // prints an arbitrary size matrix to the standard output
// void printMatrix(double **a, int rows, int cols);
// void printMatrix(double **a, int rows, int cols) {
//     int i,j;

//     for(i=0;i<rows;i++) {
//         for(j=0;j<cols;j++) {
//             printf("%.4lf ",a[i][j]);
//         }
//         printf("\n");
//     }
//     printf("\n");
// }

// // prints an arbitrary size vector to the standard output
// void printVector(double *v, int size);
// void printVector(double *v, int size) {
//     int i;

//     for(i=0;i<size;i++) {
//         printf("%.4lf ",v[i]);
//     }
//     printf("\n\n");
// }

// /*
//   Modified from Numerical Recipes in C
//   Given a matrix a[nRows][nCols], svdcmp() computes its singular value
//   decomposition, A = U * W * Vt.  A is replaced by U when svdcmp
//   returns.  The diagonal matrix W is output as a vector w[nCols].
//   V (not V transpose) is output as the matrix V[nCols][nCols].
// */
// int svdcmp(double **a, int nRows, int nCols, double *w, double **v) {
//     int flag,i,its,j,jj,k,l,nm;
//     double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

//     rv1 = (double*)malloc(sizeof(double)*nCols);
	
//     if(rv1 == NULL) {
//         printf("svdcmp(): Unable to allocate vector\n");
//         return(-1);
//     }

//     g = scale = anorm = 0.0;
//     for(i=0;i<nCols;i++) {
//         l = i+1;
//         rv1[i] = scale*g;
//         g = s = scale = 0.0;
//         if(i < nRows) {
//             for(k=i;k<nRows;k++) scale += fabs(a[k][i]);
//             if(scale) {
//                 for(k=i;k<nRows;k++) {
//                     a[k][i] /= scale;
//                     s += a[k][i] * a[k][i];
//                 }
//                 f = a[i][i];
//                 g = -SIGN(sqrt(s),f);
//                 h = f * g - s;
//                 a[i][i] = f - g;
//                 for(j=l;j<nCols;j++) {
//                     for(s=0.0,k=i;k<nRows;k++) s += a[k][i] * a[k][j];
//                     f = s / h;
//                     for(k=i;k<nRows;k++) a[k][j] += f * a[k][i];
//                 }
//                 for(k=i;k<nRows;k++) a[k][i] *= scale;
//             }
//         }
//         w[i] = scale * g;
//         g = s = scale = 0.0;
//         if(i < nRows && i != nCols-1) {
//             for(k=l;k<nCols;k++) scale += fabs(a[i][k]);
//             if(scale)  {
//                 for(k=l;k<nCols;k++) {
//                     a[i][k] /= scale;
//                     s += a[i][k] * a[i][k];
//                 }
//                 f = a[i][l];
//                 g = - SIGN(sqrt(s),f);
//                 h = f * g - s;
//                 a[i][l] = f - g;
//                 for(k=l;k<nCols;k++) rv1[k] = a[i][k] / h;
//                 for(j=l;j<nRows;j++) {
//                     for(s=0.0,k=l;k<nCols;k++) s += a[j][k] * a[i][k];
//                     for(k=l;k<nCols;k++) a[j][k] += s * rv1[k];
//                 }
//                 for(k=l;k<nCols;k++) a[i][k] *= scale;
//             }
//         }
//         anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));

//         printf(".");
//         fflush(stdout);
//     }

//     for(i=nCols-1;i>=0;i--) {
//         if(i < nCols-1) {
//             if(g) {
//                 for(j=l;j<nCols;j++)
//                     v[j][i] = (a[i][j] / a[i][l]) / g;
//                 for(j=l;j<nCols;j++) {
//                     for(s=0.0,k=l;k<nCols;k++) s += a[i][k] * v[k][j];
//                     for(k=l;k<nCols;k++) v[k][j] += s * v[k][i];
//                 }
//             }
//             for(j=l;j<nCols;j++) v[i][j] = v[j][i] = 0.0;
//         }
//         v[i][i] = 1.0;
//         g = rv1[i];
//         l = i;
//         printf(".");
//         fflush(stdout);
//     }

//     for(i=IMIN(nRows,nCols) - 1;i >= 0;i--) {
//         l = i + 1;
//         g = w[i];
//         for(j=l;j<nCols;j++) a[i][j] = 0.0;
//         if(g) {
//             g = 1.0 / g;
//             for(j=l;j<nCols;j++) {
//                 for(s=0.0,k=l;k<nRows;k++) s += a[k][i] * a[k][j];
//                 f = (s / a[i][i]) * g;
//                 for(k=i;k<nRows;k++) a[k][j] += f * a[k][i];
//             }
//             for(j=i;j<nRows;j++) a[j][i] *= g;
//         }
//         else
//             for(j=i;j<nRows;j++) a[j][i] = 0.0;
//         ++a[i][i];
//         printf(".");
//         fflush(stdout);
//     }

//     for(k=nCols-1;k>=0;k--) {
//         for(its=0;its<30;its++) {
//             flag = 1;
//             for(l=k;l>=0;l--) {
//                 nm = l-1;
//                 if((fabs(rv1[l]) + anorm) == anorm) {
//                     flag =  0;
//                     break;
//                 }
//                 if((fabs(w[nm]) + anorm) == anorm) break;
//             }
//             if(flag) {
//                 c = 0.0;
//                 s = 1.0;
//                 for(i=l;i<=k;i++) {
//                     f = s * rv1[i];
//                     rv1[i] = c * rv1[i];
//                     if((fabs(f) + anorm) == anorm) break;
//                     g = w[i];
//                     double* h1;
//                     double* f1 = &f;
//                     double* g1 = &g;
//                     svd(3, f1, g1, h1);
//                     h = *h1;
//                     w[i] = h;
//                     h = 1.0 / h;
//                     c = g * h;
//                     s = -f * h;
//                     for(j=0;j<nRows;j++) {
//                         y = a[j][nm];
//                         z = a[j][i];
//                         a[j][nm] = y * c + z * s;
//                         a[j][i] = z * c - y * s;
//                     }
//                 }
//             }
//             z = w[k];
//             if(l == k) {
//                 if(z < 0.0) {
//                     w[k] = -z;
//                     for(j=0;j<nCols;j++) v[j][k] = -v[j][k];
//                 }
//                 break;
//             }
//             if(its == 29) printf("no convergence in 30 svdcmp iterations\n");
//             x = w[l];
//             nm = k-1;
//             y = w[nm];
//             g = rv1[nm];
//             h = rv1[k];
//             f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
//             double* f1 = &f;
//             double nums = 1.0;
//             double* nums1 = &nums;
//             double* g1;
//             svd(3, f1, nums1, g1);
//             g = *g1;
//             f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g,f))) - h)) / x;
//             c = s = 1.0;
//             for(j=l;j<=nm;j++) {
//                 i = j+1;
//                 g = rv1[i];
//                 y = w[i];
//                 h = s * g;
//                 g = c * g;
//                 double* z1;
//                 f1 = &f;
//                 double* h1 = &h;
//                 svd(3, f1, h1, z1);
//                 z = *z1;
//                 rv1[j] = z;
//                 c = f/z;
//                 s = h/z;
//                 f = x * c + g * s;
//                 g = g * c - x * s;
//                 h = y * s;
//                 y *= c;
//                 for(jj=0;jj<nCols;jj++) {
//                     x = v[jj][j];
//                     z = v[jj][i];
//                     v[jj][j] = x * c + z * s;
//                     v[jj][i] = z * c - x * s;
//                 }
//                 f1 = &f;
//                 h1 = &h;
//                 svd(3, f1, h1, z1);
//                 z = *z1;
//                 w[j] = z;
//                 if(z) {
//                     z = 1.0 / z;
//                     c = f * z;
//                     s = h * z;
//                 }
//                 f = c * g + s * y;
//                 x = c * y - s * g;
//                 for(jj=0;jj < nRows;jj++) {
//                     y = a[jj][j];
//                     z = a[jj][i];
//                     a[jj][j] = y * c + z * s;
//                     a[jj][i] = z * c - y * s;
//                 }
//             }
//             rv1[l] = 0.0;
//             rv1[k] = f;
//             w[k] = x;
//         }
//         printf(".");
//         fflush(stdout);
//     }
//     printf("\n");

//     free(rv1);

//     return(0);
// }



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}


int dsvd(Matrix& a, int m, int n, Matrix& w, Matrix& v)
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }
  
    rv1 = (double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs((double)a[k][i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k][i] = (double)((double)a[k][i]/scale);
                    s += ((double)a[k][i] * (double)a[k][i]);
                }
                f = (double)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (double)(f - g);
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += (double)(f * (double)a[k][i]);
                    }
                }
                for (k = i; k < m; k++) 
                    a[k][i] = (double)((double)a[k][i]*scale);
            }
        }
        w[i][0] = (double)(scale * g);
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs((double)a[i][k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i][k] = (double)((double)a[i][k]/scale);
                    s += ((double)a[i][k] * (double)a[i][k]);
                }
                f = (double)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (double)(f - g);
                for (k = l; k < n; k++) 
                    rv1[k] = (double)a[i][k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += ((double)a[j][k] * (double)a[i][k]);
                        for (k = l; k < n; k++) 
                            a[j][k] += (double)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] = (double)((double)a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((double)w[i][0]) + fabs(rv1[i])));
    }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v[j][i] = (double)(((double)a[i][j] / (double)a[i][l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += ((double)a[i][k] * (double)v[k][j]);
                    for (k = l; k < n; k++) 
                        v[k][j] += (double)(s * (double)v[k][i]);
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
  
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = (double)w[i][0];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i][j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += ((double)a[k][i] * (double)a[k][j]);
                    f = (s / (double)a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += (double)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] = (double)((double)a[j][i]*g);
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm][0]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = (double)w[i][0];
                        h = PYTHAG(f, g);
                        w[i][0] = (double)h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (double)(y * c + z * s);
                            a[j][i] = (double)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k][0];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k][0] = (double)(-z);
                    for (j = 0; j < n; j++) 
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = (double)w[l][0];
            nm = k - 1;
            y = (double)w[nm][0];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i][0];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (double)(x * c + z * s);
                    v[jj][i] = (double)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j][0] = (double)z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (double)(y * c + z * s);
                    a[jj][i] = (double)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k][0] = (double)x;
        }
    }
	
    free((void*) rv1);
    return(1);
}



template <typename T> double sgn(T val)
{
    return (val > T(0)) - (val < T(0));
}

int squareSvd (Matrix a, int M, int N, Matrix& w, Matrix& v){
  string T,P,Db;

  double elapsedTime,elapsedTime2;
  timeval start,end,end2;

  if(M != N) {
	  cout<<"Error: Matrix must be square";
	  return 0;
  }

  
  // double **U,**V, *S,**U_t, **V_t, **A;
  double alpha, beta, gamma, c, zeta, t,s,sub_zeta, converge;

  int acum = 0;
  int temp1, temp2;
  converge = 1.0;

  // U = new double*[N];
  // V = new double*[N];
  // U_t = new double*[N];
  // V_t = new double*[N];
  // A = new double*[N];
  // S = new double[N];

  // for(int i =0; i<N; i++){
	// 	U[i] = new double[N];
	// 	V[i] = new double[N];
	// 	U_t[i] = new double[N];
	// 	V_t[i] = new double[N];
	// 	A[i] = new double[N];
  // }

	Matrix U(N,N);
	Matrix V(N,N);
	Matrix U_t(N,N);
	Matrix V_t(N,N);
	Matrix A(N,N);
	Matrix S(N,1);

  for(int i = 0; i < M; i++){
    for(int j =0; j < N; j++){
      U_t[i][j] = a[i][j];
    }
  }
 
  for(int i=0; i<M;i++){
    for(int j=0; j<N;j++){

      if(i==j){
        V_t[i][j] = 1.0;
      }
      else{
        V_t[i][j] = 0.0;
      }
    }
  }

	for(int i=0; i<M;i++){
		for(int j=0; j<N;j++){

			A[i][j] = U_t[j][i];
		}
	}


  /* SVD using Jacobi algorithm (Sequencial)*/

  //  gettimeofday(&start, NULL);

   double conv;
   while(converge > epsilon){ 		//convergence
    converge = 0.0;	
   		
    acum++;				//counter of loops

    for(int i = 1; i<M; i++){
      for(int j = 0; j<i; j++){


          alpha = 0.0;
          beta = 0.0;
          gamma = 0.0;

          for(int k = 0; k<N ; k++){
            alpha = alpha + (U_t[i][k] * U_t[i][k]);
            beta = beta + (U_t[j][k] * U_t[j][k]);
            gamma = gamma + (U_t[i][k] * U_t[j][k]);
          }

          converge = max(converge, abs(gamma)/sqrt(alpha*beta));	//compute convergence
	  								//basicaly is the angle
									//between column i and j


          zeta = (beta - alpha) / (2.0 * gamma);
          t = sgn(zeta) / (abs(zeta) + sqrt(1.0 + (zeta*zeta)));        //compute tan of angle
          c = 1.0 / (sqrt (1.0 + (t*t)));				//extract cos
          s = c*t;							//extrac sin
 

	  //Apply rotations on U and V

  	  for(int k=0; k<N; k++){
            t = U_t[i][k];
            U_t[i][k] = c*t - s*U_t[j][k];
            U_t[j][k] = s*t + c*U_t[j][k];

            t = V_t[i][k];
            V_t[i][k] = c*t - s*V_t[j][k];
            V_t[j][k] = s*t + c*V_t[j][k];

          }

      }
    }
 }


  //Create matrix S

  for(int i =0; i<M; i++){

    t=0;
    for(int j=0; j<N;j++){
      t=t + pow(U_t[i][j],2);
    }
    t = sqrt(t);

    for(int j=0; j<N;j++){
      U_t[i][j] = U_t[i][j] / t;
      if(i == j){
        S[i][0] = t;
      }
    }
  }

  // gettimeofday(&end, NULL);
 /************************************************************/

 /* Develop SVD Using OpenMP */

	// fix final result
  for(int i =0; i<M; i++){
    for(int j =0; j<N; j++){
      U[i][j] = U_t[j][i];
      V[i][j] = V_t[j][i];
    }
  }

	printMatrix(U);
	printMatrix(S);
	printMatrix(V);
  //Output time and iterations


	// cout<<"iterations: "<<acum<<endl;
	// elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
	// elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
	// cout<<"Time: "<<elapsedTime<<" ms."<<endl<<endl;




  // // Output the matrixes for debug
  // cout<<"U"<<endl<<endl;
  // for(int i =0; i<M; i++){
  //   for(int j =0; j<N; j++){

  //     cout<<U[i][j]<<"  ";
  //   }
  //   cout<<endl;
  // }

  // cout<<endl<<"V"<<endl<<endl;
  // for(int i =0; i<M; i++){
  //   for(int j =0; j<N; j++){

  //     cout<<V[i][j]<<"  ";
  //   }
  //   cout<<endl;
  // }

  // cout<<endl<<"S"<<endl<<endl;
  // for(int i =0; i<M; i++){
  //   for(int j =0; j<N; j++){

  //      if(i==j){  cout<<S[i]<<"  ";}
	
  //      else{
	//        cout<<"0.0  ";
  //      }
  //   }
  //   cout<<endl;
  // }


  //  delete [] S;
  //  for(int i = 0; i<N;i++){
	//    delete [] A[i];
	//    delete [] U[i];
	//    delete [] V[i];
	//    delete [] U_t[i];
	//    delete [] V_t[i];
  //  }

  return 0;
}

