lib_LTLIBRARIES = libgslsp.la
libgslsp_la_SOURCES = \
  spcompress.c        \
  spdgemv.c           \
  spgetset.c          \
  spmatrix.c          \
  spoper.c            \
	spswap.c

check_PROGRAMS = test
test_SOURCES = test.c test_dgemv.c test_getset.c

TESTS = $(check_PROGRAMS)

pkginclude_HEADERS = gsl_spmatrix.h

test_LDADD = libgslsp.la -lgsl -lgslcblas -lm
