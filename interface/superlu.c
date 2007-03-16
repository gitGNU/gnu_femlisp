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
 Compiling the interface file:
  gcc -c -fPIC -I/usr/include/superlu superlu.c
  ld -lsuperlu -shared superlu.o -o superlu.so
   
 Compiling for obtaining a filter:
  gcc -D__FILTER__ -g  -I/usr/include/superlu -lsuperlu -o superlu superlu.c

 For testing the direct call:
  gcc -D__TESTING__  -I/usr/include/superlu -lsuperlu superlu.c

 For testing the filter:
  gcc -D__TEST_FILTER__ -g superlu.c
************************************************************************/

#ifndef  __TEST_FILTER__
#include "dsp_defs.h"

int c_superlu (int m, int n, int nnz, int *Ap, int *Ai, double *Ax,
               int nrhs, double *rhs, double *sol, int orientation)
{
    SuperMatrix A, B, L, U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int i, info;
    superlu_options_t options;
    SuperLUStat_t stat;

    if ((orientation!=0)&&(orientation!=1)) return -1;
    
    /* copy rhs to solution */
    if (sol!=rhs)
        for (i=0; i<m*nrhs; i++) sol[i] = rhs[i];
    
    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix(&A, m, n, nnz, Ax, Ai, Ap,
                           (orientation==0)?SLU_NC:SLU_NR,
                           SLU_D, SLU_GE);
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
#else

#include <stdio.h>
#include <stdlib.h>

int call_superlu (int m, int n, int nnz, int *Ap, int *Ai, double *Ax,
                  int nrhs, double *rhs, double *sol, int orientation)
{
    FILE *file;
    char command[255];

    // write lgs
    file=fopen("superlu.lgs","w");
    if (file==0) return -1;
    
    /* write matrix and rhs to file */
    fwrite (Ap, sizeof(int), n+1, file);
    fwrite (Ai, sizeof(int), nnz, file);
    fwrite (Ax, sizeof(double), nnz, file);
    fwrite (rhs, sizeof(double), m*nrhs, file);
    fclose(file);

    // call solve
    sprintf(command, "cat superlu.lgs | superlu %d %d %d %d %d | cat > superlu.solution",
            m, n, nnz, nrhs, orientation);
    system(command);

    // read solution
    file=fopen("superlu.solution","r");
    if (file==0) return -2;
    fread (sol, sizeof(double), m*nrhs, file);
    fclose(file);

    return 0;
}
#endif

#ifdef __FILTER__

#include <stdio.h>

int main(int argc, char *argv[])
{
    int nrhs, m, n, nnz, orientation;
    int *asub, *xa;
    double *a, *rhs, *sol;
	
    if (argc!=5) return -1;
    if (sscanf (argv[1], "%d", &m)!=1) return -2;
    if (sscanf (argv[2], "%d", &n)!=1) return -3;
    if (sscanf (argv[3], "%d", &nnz)!=1) return -4;
    if (sscanf (argv[4], "%d", &nrhs)!=1) return -5;
    if (sscanf (argv[5], "%d", &orientation)!=1) return -6;

    if (!(xa = malloc ((n+1)*sizeof(int)))) return -6;
    if (!(asub = malloc (nnz*sizeof(int)))) return -7;
    if (!(a = malloc (nnz*sizeof(double)))) return -8;
    if (!(rhs = malloc (m*nrhs*sizeof(double)))) return -9;

    /* read matrix and rhs from stdin */
    fread (xa, sizeof(int), n+1, stdin);
    fread (asub, sizeof(int), nnz, stdin);
    fread (a, sizeof(double), nnz, stdin);
    fread (rhs, sizeof(double), m*nrhs, stdin);
	
    sol = rhs;  /* rhs is overwritten */

    if (c_superlu (m, n, nnz, xa, asub, a, nrhs, rhs, sol, orientation) != 0)
        return -10;
    
    /* write out solution */
    fwrite (rhs, sizeof(double), m*nrhs, stdout);
    
    return 0;
}

#else

void fill_vector (double *x, int n, double s)
{
    int i;
    for (i=0; i<n; i++) x[i] = s;
}

void print_vector (double *x, int n)
{
    int i;
    for (i=0; i<n; i++) printf("%lf\n", x[i]);
    printf("\n");
}

/* Test example:
   
   s=19; u=21; p=16; e=5; r=18; l= 12;
   A=[s,0,u,u,0; l,u,0,0,0; 0,l,p,0,0; 0,0,0,e,u; l,l,0,0,r];
   A\[1,1,1,1,1]'
   A'\[1,1,1,1,1]'

   Should give:
   
-0.031250
0.065476
0.013393
0.062500
0.032738

0.033298
0.045232
0.018797
0.060150
-0.014620
   
*/

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
    sol = rhs;  /* rhs is overwritten */

    fill_vector (rhs, m*nrhs, 1.0);
#ifdef __TEST_FILTER__
    call_superlu (m, n, nnz, xa, asub, a, nrhs, rhs, sol, 0);
#else
    c_superlu (m, n, nnz, xa, asub, a, nrhs, rhs, sol, 0);
#endif
    print_vector (sol,m*nrhs);
    
    fill_vector (rhs, m*nrhs, 1.0);
#ifdef __TEST_FILTER__
    call_superlu (m, n, nnz, xa, asub, a, nrhs, rhs, sol, 1);
#else
    c_superlu (m, n, nnz, xa, asub, a, nrhs, rhs, sol, 1);
#endif
    print_vector (sol, m*nrhs);
    
    return 0;
}
#endif
