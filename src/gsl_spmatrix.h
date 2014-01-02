/* gsl_spmatrix.h
 * 
 * Copyright (C) 2012 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SPMATRIX_H__
#define __GSL_SPMATRIX_H__

#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct
{
  size_t size1; /* number of rows */
  size_t size2; /* number of columns */

  size_t *i;    /* row indices of size nzmax */
  size_t *j;    /* column indices of size nzmax */
  double *data; /* matrix elements of size nzmax */

  /*
   * p contains the column or row pointers.
   *
   * comp. col: p[j] = index in data of first non-zero element in column j
   * comp. row: p[i] = index in data of first non-zero element in row i
   */
  size_t *p;

  size_t nzmax; /* maximum number of matrix elements */
  size_t nz;    /* number of non-zero values in matrix */

  size_t flags;
} gsl_spmatrix;

#define GSL_SPMATRIX_TRIPLET      (1 << 0)
#define GSL_SPMATRIX_COMPCOL      (1 << 1)
#define GSL_SPMATRIX_COMPROW      (1 << 1)

/* default value of nzmax */
#define GSL_SPMATRIX_NZMAX        10000

#define GSLSP_ISTRIPLET(m)        ((m)->flags & GSL_SPMATRIX_TRIPLET)
#define GSLSP_ISCOMPCOL(m)        ((m)->flags & GSL_SPMATRIX_COMPCOL)

/*
 * Prototypes
 */

gsl_spmatrix *gsl_spmatrix_alloc(const size_t n1, const size_t n2);
gsl_spmatrix *gsl_spmatrix_alloc_nzmax(const size_t n1, const size_t n2,
                                       const size_t nzmax);
void gsl_spmatrix_free(gsl_spmatrix *m);
int gsl_spmatrix_realloc(const size_t nzmax, gsl_spmatrix *m);
int gsl_spmatrix_reset(gsl_spmatrix *m);
size_t gsl_spmatrix_nnz(const gsl_spmatrix *m);

double gsl_spmatrix_get(const gsl_spmatrix *m, const size_t i,
                        const size_t j);
int gsl_spmatrix_set(gsl_spmatrix *m, const size_t i, const size_t j,
                     const double x);

gsl_spmatrix *gsl_spmatrix_compcol(const gsl_spmatrix *T);

/* spoper.c */
int gsl_spmatrix_scale(gsl_spmatrix *m, const double x);
int gsl_spmatrix_minmax(const gsl_spmatrix *m, double *min_out,
                        double *max_out);
int gsl_spmatrix_d2sp(gsl_spmatrix *S, const gsl_matrix *A);
int gsl_spmatrix_sp2d(gsl_matrix *A, const gsl_spmatrix *S);

/* spblas */
int gsl_spblas_dgemv(const double alpha, const gsl_spmatrix *A,
                     const gsl_vector *x, const double beta, gsl_vector *y);

__END_DECLS

#endif /* __GSL_SPMATRIX_H__ */
