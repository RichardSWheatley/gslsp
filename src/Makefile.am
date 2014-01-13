lib_LTLIBRARIES = libgslsp.la
libgslsp_la_SOURCES = \
  spcompress.c        \
	spcopy.c            \
  spdgemv.c           \
  spdgemm.c           \
  spgetset.c          \
  spmatrix.c          \
  spoper.c            \
	spprop.c            \
	spswap.c

check_PROGRAMS = test
test_SOURCES = test.c

TESTS = $(check_PROGRAMS)

pkginclude_HEADERS = gsl_spmatrix.h

test_LDADD = libgslsp.la -lgsl -lgslcblas -lm
