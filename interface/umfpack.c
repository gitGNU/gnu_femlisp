/***********************************************************************
  
 umfpack.c - Femlisp interface to UMFPACK

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
  gcc -I/usr/include/suitesparse/ -Wall -fPIC -c umfpack.c
  ld -lumfpack -lamd -shared umfpack.o -o umfpack.so
 For testing purposes:
  gcc -D__TESTING__ -I/usr/include/suitesparse/ -lumfpack -lamd umfpack.c
************************************************************************/

#include <stdlib.h>
#include "umfpack.h"

int c_umfpack (int m, int n, int nnz, int *Ap, int *Ai, double *Ax,
               int nrhs, double *B, double *X, int orientation)
{
    int i, status;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];
    void *Symbolic = 0, *Numeric = 0;

    if (X==B)
    {
        /* UMFPACK can't handle this case, so we fix it by copying B */
        if (!(B = malloc(nrhs*m*sizeof(double))))
            return UMFPACK_ERROR_out_of_memory;
        for (i=0; i<nrhs*m; i++) B[i] = X[i];
    }
    
    umfpack_di_defaults (Control);
    
    status = umfpack_di_symbolic (m, n, Ap, Ai, Ax, &Symbolic, Control, Info);
    if (status!=UMFPACK_OK) goto cleanup_and_return;
    status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
    if (status!=UMFPACK_OK) goto cleanup_and_return;
    for (i=0; i<nrhs; i++)
    {
        status = umfpack_di_solve ((orientation==0)?UMFPACK_A:UMFPACK_At,
                                   Ap, Ai, Ax, X+i*n, B+i*m, Numeric, Control, Info);
        if (status!=UMFPACK_OK) goto cleanup_and_return;
    }
    
 cleanup_and_return:
    if (Symbolic) umfpack_di_free_symbolic (&Symbolic);
    if (Numeric) umfpack_di_free_numeric (&Numeric);
    
    return status;
}

#ifdef __TESTING__
#include <stdio.h>

void fill (double *x, int n, double s)
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
    if (!(rhs = malloc(m*nrhs*sizeof(double)))) return -1;

    sol = rhs;  /* rhs is overwritten */
    
    fill(rhs, m*nrhs, 1.0);
    c_umfpack (m, n, nnz, xa, asub, a, nrhs, rhs, sol, 0);
    print_vector (sol,m*nrhs);

    fill(rhs, m*nrhs, 1.0);
    c_umfpack (m, n, nnz, xa, asub, a, nrhs, rhs, sol, 1);
    print_vector (sol,m*nrhs);

    return 0;
}
#endif
