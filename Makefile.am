
bin_PROGRAMS = wgms3d

wgms3d_SOURCES = \
	main.cc \
	wgms3d.cc \
	wgms3d.h \
	diffops.cc \
	diffops.h \
	stencil.h \
	sparse.h \
	mgp.cc \
	mgp.h \
	geometry.h \
	bezier.cc \
	bezier.h \
	fortran_interface.h \
	standarddiffop_2nd.h \
	standarddiffop_handmade.h \
	pml.h

INCLUDES = ${SUPERLU_INCLUDES} ${ARPACK_INCLUDES}

wgms3d_LDADD = ${SUPERLU_LIBS} ${ARPACK_LIBS}

wgms3d_CXXFLAGS = -Wall -std=c++0x

DIST_SUBDIRS = tests matlab

# The following is a Maple script used to generate the
# homogeneous-space differential operators (standarddiffop_2nd.h).
EXTRA_DIST = \
	standarddiffop.mpl \
	NOTES
