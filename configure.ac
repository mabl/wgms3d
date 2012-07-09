AC_INIT
AM_INIT_AUTOMAKE(wgms3d, 1.2)
AC_CONFIG_SRCDIR([wgms3d.cc])
AC_CONFIG_HEADER(config.h)

AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_WRAPPERS

AC_PROG_CXX

sinclude(acx_blas.m4)
sinclude(acx_lapack.m4)

ACX_BLAS([], [AC_MSG_ERROR([BLAS was not found!])])
ACX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found!])])

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

##############################################################################
# Check to see if calling Fortran functions (in particular, the BLAS
# complex dot product functions) works from C:
# (adapted from mpb-1.4.2)

AC_MSG_CHECKING([whether calling BLAS zdotc works])
AC_F77_FUNC(zdotc)
AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
typedef struct { double re, im; } cmplx;
#ifdef __cplusplus
extern "C"
#endif
extern cmplx $zdotc(int *, cmplx *, int *, cmplx *, int *);
int main(void) {
     cmplx x[2] = { {1,0}, {0,1} }, y[2] = { {1,0}, {1,0} };
     int n = 2, inc1 = 1, inc2 = 1;
     cmplx result = { 0, 0 };
     result = $zdotc(&n, x, &inc1, y, &inc2);
     if (result.re != 1 || result.im != -1) return EXIT_FAILURE;
     else return EXIT_SUCCESS;
}
], ok="yes", ok="no", ok="no")
AC_MSG_RESULT($ok)
if test "$ok" = "no"; then
	AC_MSG_ERROR(Your BLAS doesn't work.)
fi

##############################################################################
# Check for SuperLU
#

AC_ARG_WITH(superlu,
	    [  --with-superlu=DIR	  compile against SuperLU installed in DIR],
	    [
	        SUPERLU_LIBS="${withval}/lib/libsuperlu.so.4.0.0"
		SUPERLU_INCLUDES="-I${withval}/include"
	    ])
AC_SUBST(SUPERLU_LIBS)
AC_SUBST(SUPERLU_INCLUDES)

save_LIBS="$LIBS";
save_CFLAGS="$CFLAGS";
LIBS="$SUPERLU_LIBS $LIBS"
CFLAGS="$SUPERLU_INCLUDES $CFLAGS"
AC_MSG_CHECKING([whether linking to SuperLU works])
AC_TRY_LINK([
	#include <slu_ddefs.h>
	#include <slu_util.h>
	],[
    superlu_options_t *a;
    SuperMatrix *AC;
    int d, e, *f;
    void *g;
    int h, *i, *j;
    SuperMatrix *k, *l;
    SuperLUStat_t *m;
    int *n;
    dgstrf(a,AC,d,e,f,g,h,i,j,k,l,m,n);
    F77_FUNC(dcopy,DCOPY)();
], ok="yes", ok="no", ok="no")
LIBS="$save_LIBS"
CFLAGS="$save_CFLAGS"
AC_MSG_RESULT($ok)
if test "$ok" = "no"; then
	AC_MSG_ERROR(Please try --with-superlu={prefix}.)
fi

##############################################################################
# Check for ARPACK
#

AC_ARG_WITH(arpack,
	    [  --with-arpack=DIR	  compile against ARPACK installed in DIR],
	    [
	        ARPACK_LIBS="${withval}/lib/libarpack.so.96.0.3"
		ARPACK_INCLUDES="-I${withval}/include"
	    ])
AC_SUBST(ARPACK_LIBS)
AC_SUBST(ARPACK_INCLUDES)

save_LIBS="$LIBS";
save_CFLAGS="$CFLAGS";
LIBS="$ARPACK_LIBS $LIBS"
CFLAGS="$ARPACK_INCLUDES $CFLAGS"
AC_MSG_CHECKING([whether linking to ARPACK works])
AC_F77_FUNC(dnaupd)
AC_TRY_LINK([],[
    $dnaupd();
    F77_FUNC(dcopy,DCOPY)();
], ok="yes", ok="no", ok="no")
LIBS="$save_LIBS"
CFLAGS="$save_CFLAGS"
AC_MSG_RESULT($ok)
if test "$ok" = "no"; then
	AC_MSG_ERROR(Please try --with-arpack={prefix}.)
fi

AC_CONFIG_FILES([Makefile
		 matlab/Makefile
		 tests/Makefile
		 tests/all_field_components/Makefile
		 tests/complex_modes/Makefile
		 tests/fiber_convergence/Makefile
		 tests/lossy_materials/Makefile
		 tests/semivectorial/Makefile
		 tests/disk_resonator/Makefile
		 tests/silicon_strip_waveguide/Makefile
		 tests/tm2te_leakage/Makefile])
AC_OUTPUT
