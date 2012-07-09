
/*
    wgms3d - a full-vectorial finite-difference mode solver.

    Copyright (C) 2005-2012  Michael Krause <m.krause@tu-harburg.de>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _FORTRAN_INTERFACE_H
#define _FORTRAN_INTERFACE_H

#include "config.h"

#include <slu_ddefs.h>
#include <slu_util.h>

#include <complex>
using std::complex;

/* The following works on Linux, and on HPUX in 32-bit and 64-bit modes */
typedef int logical;

extern "C" {

    /* ARPACK functions */
    extern void
    F77_FUNC(znaupd,ZNAUPD) (int *ido, char *bmat, int *n, char *which, int *nev,
			     double *tol, std::complex<double> *resid, int *ncv, std::complex<double> *v,
			     int *ldv, int *iparam, int *ipntr, std::complex<double> *workd,
			     complex<double> *workl, int *lworkl, double *rwork, int *info);
    extern void
    F77_FUNC(dnaupd,DNAUPD) (int *ido, char *bmat, int *n, char *which, int *nev,
			     double *tol, double *resid, int *ncv, double *v,
			     int *ldv, int *iparam, int *ipntr, double *workd,
			     double *workl, int *lworkl, int *info);
    extern void
    F77_FUNC(zneupd,ZNEUPD) (logical *rvec, char *howmny, logical *select, complex<double> *d,
			     complex<double> *z, int *ldz, complex<double> *sigma, complex<double> *workev,
			     char *bmat, int *n, char *which, int *nev, double *tol,
			     complex<double> *resid, int *ncv, complex<double> *v, int *ldv,
			     int *iparam, int *ipntr, complex<double> *workd, complex<double> *workl,
			     int *lworkl, double *rwork, int *info);
    extern void
    F77_FUNC(dneupd,DNEUPD) (logical *rvec, char *howmny, logical *select, double *dr, double *di,
			     double *z, int *ldz, double *sigmar, double *sigmai, double *workev,
			     char *bmat, int *n, char *which, int *nev, double *tol, double *resid,
			     int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd,
			     double *workl, int *lworkl, int *info);

    /* BLAS/LAPACK functions */
    extern void F77_FUNC(zgels,ZGELS) (char *TRANS, int *M, int *N, int *NRHS, complex<double> *A, int *LDA,
			complex<double> *B, int *LDB, complex<double> *WORK, int *LWORK, int *INFO);
    extern void F77_FUNC(dcopy,DCOPY) (int *n, double *x, int *incx, double *y, int *incy);
    extern void F77_FUNC(zcopy,ZCOPY) (int *n, complex<double> *x, int *incx, complex<double> *y, int *incy);
    extern void F77_FUNC(dscal,DSCAL) (int *n, double *alpha, double *x, int *incx);
    extern void F77_FUNC(zscal,ZSCAL) (int *n, complex<double> *alpha, complex<double> *x, int *incx);
    extern void F77_FUNC(daxpy,DAXPY) (int *n, double *alpha, double *x,
			int *incx, double *y, int *incy);
    extern void F77_FUNC(zaxpy,ZAXPY) (int *n, complex<double> *alpha, complex<double> *x,
			int *incx, complex<double> *y, int *incy);
    extern void F77_FUNC(dgemm,DGEMM) (char *TRANSA, char *TRANSB, int *M, int *N, int *K,
			double *ALPHA, double *A, int *LDA, double *B, int *LDB,
			double *BETA, double *C, int *LDC);
    extern void F77_FUNC(zgemm,ZGEMM) (char *TRANSA, char *TRANSB, int *M, int *N, int *K,
			complex<double> *ALPHA, complex<double> *A, int *LDA,
			complex<double> *B, int *LDB,
			complex<double> *BETA, complex<double> *C, int *LDC);

    /* SuperLU functions */
    /* We want to use both the real and complex SuperLU routines, but we
       can't include slu_ddefs.h and slu_zdefs.h at the same time, since
       they contain colliding definitions.

       Therefore we copy here the definitions from slu_zdefs.h...

       This is for SuperLU Version 4.0. */

    extern void    zgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
			   SuperMatrix *, SuperLUStat_t*, int *);
    
    extern void    zgstrf (superlu_options_t*, SuperMatrix*,
			   int, int, int*, void *, int, int *, int *, 
			   SuperMatrix *, SuperMatrix *, SuperLUStat_t*, int *);

    extern void
    zCreate_CompCol_Matrix(SuperMatrix *, int, int, int, complex<double> *,
			   int *, int *, Stype_t, Dtype_t, Mtype_t);

    extern int     
    zQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
}

namespace {

    template <class T> T sq (T x) {
	return x * x;
    }

    template <class T> T re (T &x) {
	return x;
    }
    
    template <class T> T re (complex<T> &x) {
	return x.real();
    }
    
    template <class T> T im (T &x) {
	return 0.0;
    }

    template <class T> T im (complex<T> &x) {
	return x.imag();
    }
    
    void COPY (int n, double *x, int incx, double *y, int incy) {
	F77_FUNC(dcopy,DCOPY) (&n, x, &incx, y, &incy); }
    void COPY (int n, complex<double> *x, int incx, complex<double> *y, int incy) {
	F77_FUNC(zcopy,ZCOPY) (&n, x, &incx, y, &incy); }
    void SCAL (int n, double alpha, double *x, int incx) {
	F77_FUNC(dscal,DSCAL) (&n, &alpha, x, &incx); }
    void SCAL (int n, complex<double> alpha, complex<double> *x, int incx) {
	F77_FUNC(zscal,ZSCAL) (&n, &alpha, x, &incx); }
    void GELS (char *TRANS, int *M, int *N, int *NRHS, complex<double> *A, int *LDA,
	       complex<double> *B, int *LDB, complex<double> *WORK, int *LWORK, int *INFO) {
	F77_FUNC(zgels,ZGELS) (TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO); }
    void AXPY (int n, double alpha, double *x,
		      int incx, double *y, int incy) {
	F77_FUNC(daxpy,DAXPY) (&n, &alpha, x, &incx, y, &incy); }
    void AXPY (int n, complex<double> alpha, complex<double> *x,
		      int incx, complex<double> *y, int incy) {
	F77_FUNC(zaxpy,ZAXPY) (&n, &alpha, x, &incx, y, &incy); }
    void GEMM (char *TRANSA, char *TRANSB, int *M, int *N, int *K,
		      double *ALPHA, double *A, int *LDA, double *B, int *LDB,
		      double *BETA, double *C, int *LDC) {
	F77_FUNC(dgemm,DGEMM) (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC); }
    void GEMM (char *TRANSA, char *TRANSB, int *M, int *N, int *K,
		      complex<double> *ALPHA, complex<double> *A, int *LDA,
		      complex<double> *B, int *LDB,
		      complex<double> *BETA, complex<double> *C, int *LDC) {
	F77_FUNC(zgemm,ZGEMM) (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC); }

    void
	NAUPD (int *ido, char *bmat, int *n, char *which, int *nev,
	       double *tol, complex<double> *resid, int *ncv, complex<double> *v,
	       int *ldv, int *iparam, int *ipntr, complex<double> *workd,
	       complex<double> *workl, int *lworkl, double *rwork, int *info)
	{
	    F77_FUNC(znaupd,ZNAUPD) (ido, bmat, n, which, nev, tol, resid, ncv, v,
			 ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
	}
    
    void
	NAUPD (int *ido, char *bmat, int *n, char *which, int *nev,
	       double *tol, double *resid, int *ncv, double *v,
	       int *ldv, int *iparam, int *ipntr, double *workd,
	       double *workl, int *lworkl, double *rwork, int *info)
	{
	    F77_FUNC(dnaupd,DNAUPD) (ido, bmat, n, which, nev, tol, resid, ncv, v,
			 ldv, iparam, ipntr, workd, workl, lworkl, info);
	}
    
    void
	NEUPD (logical *rvec, char *howmny, logical *select, complex<double> *d,
	       complex<double> *z, int *ldz, complex<double> sigma, complex<double> *workev,
	       char *bmat, int *n, char *which, int *nev, double *tol,
	       complex<double> *resid, int *ncv, complex<double> *v, int *ldv,
	       int *iparam, int *ipntr, complex<double> *workd, complex<double> *workl,
	       int *lworkl, double *rwork, int *info)
	{
	    F77_FUNC(zneupd,ZNEUPD) (rvec, howmny, select, d, z, ldz, &sigma,
			 workev, bmat, n, which, nev, tol, resid, ncv,
			 v, ldv, iparam, ipntr, workd, workl, lworkl,
			 rwork, info);
	}
    
    void
	NEUPD (logical *rvec, char *howmny, logical *select, complex<double> *d,
	       double *z, int *ldz, double sigma, double *workev,
	       char *bmat, int *n, char *which, int *nev, double *tol,
	       double *resid, int *ncv, double *v, int *ldv,
	       int *iparam, int *ipntr, double *workd, double *workl,
	       int *lworkl, double *rwork, int *info)
	{
	    int i;
	    int nevcopy = *nev;

	    double *dr = new double[nevcopy+1];
	    double *di = new double[nevcopy+1];
	    
	    double sigmar = sigma;
	    double sigmai = 0.0;
	    
	    F77_FUNC(dneupd,DNEUPD) (rvec, howmny, select, dr, di, z, ldz, &sigmar, &sigmai,
			 workev, bmat, n, which, nev, tol, resid, ncv,
			 v, ldv, iparam, ipntr, workd, workl, lworkl, info);
	    
	    for(i = 0; i < nevcopy+1; i++) {
		d[i] = complex<double>(dr[i],di[i]);
	    }
	    
	    delete[] dr;
	    delete[] di;
	}

#ifndef USE_PARDISO
    void
	Create_CompCol_Matrix (SuperMatrix *a, int b, int c, int d, double *e,
			       int *f, int *g, Stype_t h, Mtype_t j)
	{
	    dCreate_CompCol_Matrix (a, b, c, d, e, f, g, h, SLU_D, j);
	}
    
    void
	Create_CompCol_Matrix (SuperMatrix *a, int b, int c, int d, complex<double> *e,
			       int *f, int *g, Stype_t h, Mtype_t j)
	{
	    zCreate_CompCol_Matrix (a, b, c, d, reinterpret_cast <complex<double>*> (e), f, g, h, SLU_Z, j);
	}

    Dtype_t
	get_slu_type (double *)
	{
	    return SLU_D;
	}
    Dtype_t
	get_slu_type (complex<double> *)
	{
	    return SLU_Z;
	}
#endif

    int is_complex_calculation (const double *) {
	return 0;
    }

    int is_complex_calculation (const complex<double> *) {
	return 1;
    }
}

template <class T>
static void
dump_fortran_matrix (T *A,
		     int lda,
		     int nrows,
		     int ncols)
{
    int i, j;

    for(i = 0; i < nrows; i++) {
	for(j = 0; j < ncols; j++) {
	    //	    printf("% 2.1e+i% 2.1e ", A[i+j*lda].real(), A[i+j*lda].imag());
	    if(MABS(A[i+j*lda]) < 1*1e-12) {
		printf("          ");
	    } else {
		printf("% 2.2e ", A[i+j*lda]);
	    }
	}
	printf("\n");
    }
}

template <class T>
static void
dump_fortran_matrix (complex<T> *A,
		     int lda,
		     int nrows,
		     int ncols)
{
    int i, j;

    for(i = 0; i < nrows; i++) {
	for(j = 0; j < ncols; j++) {
	    printf("% 2.1e+i% 2.1e ", A[i+j*lda].real(), A[i+j*lda].imag());
	}
	printf("\n");
    }
}

template <class T>
static void
dump_fortran_matrix (T *A,
		     int nrows,
		     int ncols)
{
    dump_fortran_matrix(A, nrows, nrows, ncols);
}

template <class T>
static void
matrix_copy (T *dst,
	     int ldst,
	     T *src,
	     int lsrc,
	     int nrows,
	     int ncols)
{
    int i, j;

    for(i = 0; i < nrows; i++) {
	for(j = 0; j < ncols; j++) {
	    dst[i + j*ldst] = src[i + j*lsrc];
	}
    }
}

#endif // _FORTRAN_INTERFACE_H

