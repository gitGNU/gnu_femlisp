/***********************************************************************
  
 superlu.c - Femlisp interface to SuperLU

 Copyright (C) 2004 Nicolas Neuss, University of Heidelberg.
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:
 
 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
 WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
 NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF HEIDELBERG OR OTHER
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

************************************************************************/

/************************************************************************
 Compiling:
  gcc -c -fPIC -I/usr/include/superlu superlu.c
  ld -lsuperlu -shared superlu.o -o superlu.so
   
 Compiling for obtaining a filter:
  gcc -D__FILTER__ -g  -I/usr/include/superlu -lsuperlu -o superlu superlu.c

 For testing purposes:
  gcc -D__TESTING__  -I/usr/include/superlu -lsuperlu superlu.c
************************************************************************/

#include "dsp_defs.h"

int c_superlu (int m, int n, int nnz, int *Ap, int *Ai, double *Ax,
			   int nrhs, double *rhs, double *sol)
{
    SuperMatrix A, B, L, U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int i, info;
    superlu_options_t options;
    SuperLUStat_t stat;

	/* copy rhs to solution */
	if (sol!=rhs)
		for (i=0; i<m*nrhs; i++) sol[i] = rhs[i];
	
    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix(&A, m, n, nnz, Ax, Ai, Ap, SLU_NC, SLU_D, SLU_GE);
    /* Create right-hand side matrix B (same as X). */
    dCreate_Dense_Matrix(&B, m, nrhs, sol, m, SLU_DN, SLU_D, SLU_GE);
	
    if (!(perm_r = intMalloc(m))) return -1;
    if (!(perm_c = intMalloc(n))) {SUPERLU_FREE (perm_r); return -1;}
	
    set_default_options(&options);  /* options.ColPerm = NATURAL;*/
    StatInit(&stat);
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    
    /* De-allocate storage */
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
	
	return info;
}

#ifdef __FILTER__
#include <stdio.h>

int main(int argc, char *argv[])
{
    int nrhs, m, n, nnz;
    int *asub, *xa;
    double *a, *rhs, *sol;
	
	if (argc!=5) return -1;
	if (sscanf (argv[1], "%d", &m)!=1) return -2;
	if (sscanf (argv[2], "%d", &n)!=1) return -3;
	if (sscanf (argv[3], "%d", &nnz)!=1) return -4;
	if (sscanf (argv[4], "%d", &nrhs)!=1) return -5;

    if (!(xa = malloc ((n+1)*sizeof(int)))) return -6;
    if (!(asub = malloc (nnz*sizeof(int)))) return -7;
    if (!(a = malloc (nnz*sizeof(double)))) return -8;
    if (!(rhs = malloc (m*nrhs*sizeof(double)))) return -9;

	/* read matrix and rhs from stdin */
	fread (xa, sizeof(int), n+1, stdin);
	fread (asub, sizeof(int), nnz, stdin);
	fread (xa, sizeof(double), nnz, stdin);
	fread (rhs, sizeof(double), m*nrhs, stdin);
	
	sol = rhs;  /* rhs is overwritten */
	
	c_superlu (m, n, nnz, xa, asub, a, nrhs, rhs, sol);

	/* write out solution */
	fwrite (rhs, sizeof(double), m*nrhs, stdin);
	
	return 0;
}
#endif

#ifdef __TESTING__
#include <stdio.h>

int main(int argc, char *argv[])
{
    double   *a, *rhs, *sol;
    double   s, u, p, e, r, l;
    int      *asub, *xa;
    int      nrhs, i, m, n, nnz;

    m = n = 5;
    nnz = 12;
    if (!(a = malloc(nnz*sizeof(double)))) return -1;
    if (!(asub = malloc(nnz*sizeof(int)))) return -1;
    if (!(xa = malloc((n+1)*sizeof(int)))) return -1;
    s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
    a[0] = s; a[1] = l; a[2] = l; a[3] = u; a[4] = l; a[5] = l;
    a[6] = u; a[7] = p; a[8] = u; a[9] = e; a[10]= u; a[11]= r;
    asub[0] = 0; asub[1] = 1; asub[2] = 4; asub[3] = 1;
    asub[4] = 2; asub[5] = 4; asub[6] = 0; asub[7] = 2;
    asub[8] = 0; asub[9] = 3; asub[10]= 3; asub[11]= 4;
    xa[0] = 0; xa[1] = 3; xa[2] = 6; xa[3] = 8; xa[4] = 10; xa[5] = 12;

    nrhs = 1;
    if (!(rhs = malloc(m*nrhs*sizeof(double))))	return -1;
    for (i = 0; i < m*nrhs; ++i) rhs[i] = 1.0;
	sol = rhs;  /* rhs is overwritten */
	
	c_superlu (m, n, nnz, xa, asub, a, nrhs, rhs, sol);

	for (i=0; i<m; i++)	printf("%lf\n", sol[i]);
	return 0;
}
#endif
