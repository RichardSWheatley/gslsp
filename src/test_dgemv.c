/* test_dgemv.c
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
 *
 * This program tests the gsl_spmatrix_set() and gsl_spmatrix_get()
 * routines by generating random sparse matrices, calling the two
 * functions and comparing results
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_blas.h>

#include "gsl_spmatrix.h"
#include "test.h"

int
test_dgemv_proc(const double alpha, const double beta)
{
  int s = 0;
  size_t N_max = 100;
  size_t N, M;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_spmatrix *St = gsl_spmatrix_alloc(1, 1);
  gsl_spmatrix *Sc;

  for (M = 1; M <= N_max; ++M)
    {
      gsl_vector *y_gsl = gsl_vector_alloc(M);
      gsl_vector *y_sp = gsl_vector_alloc(M);

      for (N = 1; N <= N_max; ++N)
        {
          gsl_matrix *A = gsl_matrix_calloc(M, N);
          gsl_vector *x = gsl_vector_alloc(N);

          /* create random matrix */
          test_random_sparse_matrix(A, r, -10.0, 10.0);

          /* create random vector */
          test_random_vector(x, r, -10.0, 10.0);

          gsl_spmatrix_d2sp(St, A);
          Sc = gsl_spmatrix_compcol(St);

          /* compute y = A x with gsl */
          gsl_vector_set_zero(y_gsl);
          gsl_blas_dgemv(CblasNoTrans, alpha, A, x, beta, y_gsl);

          /* compute y = A x with spblas */
          gsl_vector_set_zero(y_sp);
          gsl_spblas_dgemv(alpha, St, x, beta, y_sp);
          test_vectors(y_sp, y_gsl, 1.0e-9, "test_dgemv: triplet format");

          gsl_vector_set_zero(y_sp);
          gsl_spblas_dgemv(alpha, Sc, x, beta, y_sp);
          test_vectors(y_sp, y_gsl, 1.0e-9, "test_dgemv: compressed column format");

          gsl_matrix_free(A);
          gsl_vector_free(x);
          gsl_spmatrix_free(Sc);
        }

      gsl_vector_free(y_gsl);
      gsl_vector_free(y_sp);
    }

  gsl_rng_free(r);
  gsl_spmatrix_free(St);

  return s;
} /* test_dgemv_proc() */

int
test_dgemv(void)
{
  int s = 0;

  s += test_dgemv_proc(1.0, 0.0);
  s += test_dgemv_proc(2.4, -0.5);
  s += test_dgemv_proc(0.1, 10.0);

  return s;
} /* test_dgemv() */
