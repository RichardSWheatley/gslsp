/* test_dgemv.c
 * 
 * Copyright (C) 2012-2014 Patrick Alken
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
  gsl_matrix *A = gsl_matrix_alloc(N_max, N_max);
  gsl_vector *x = gsl_vector_alloc(N_max);
  gsl_vector *y1 = gsl_vector_alloc(N_max);
  gsl_vector *y2 = gsl_vector_alloc(N_max);
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_spmatrix *St = gsl_spmatrix_alloc(1, 1);
  gsl_spmatrix *Sc;
  size_t N, M;

  for (M = 1; M <= N_max; ++M)
    {
      gsl_vector_view y_gsl = gsl_vector_subvector(y1, 0, M);
      gsl_vector_view y_sp = gsl_vector_subvector(y2, 0, M);

      for (N = 1; N <= N_max; ++N)
        {
          gsl_matrix_view Av = gsl_matrix_submatrix(A, 0, 0, M, N);
          gsl_vector_view xv = gsl_vector_subvector(x, 0, N);

          /* create random matrix */
          test_random_sparse_matrix(&Av.matrix, r, -10.0, 10.0);

          /* create random vector */
          test_random_vector(&xv.vector, r, -10.0, 10.0);

          gsl_spmatrix_d2sp(St, &Av.matrix);
          Sc = gsl_spmatrix_compcol(St);

          /* compute y = A x with gsl */
          gsl_vector_set_zero(&y_gsl.vector);
          gsl_blas_dgemv(CblasNoTrans, alpha, &Av.matrix, &xv.vector, beta, &y_gsl.vector);

          /* compute y = A x with spblas */
          gsl_vector_set_zero(&y_sp.vector);
          gsl_spblas_dgemv(alpha, St, &xv.vector, beta, &y_sp.vector);
          test_vectors(&y_sp.vector, &y_gsl.vector, 1.0e-9,
                       "test_dgemv: triplet format");

          gsl_vector_set_zero(&y_sp.vector);
          gsl_spblas_dgemv(alpha, Sc, &xv.vector, beta, &y_sp.vector);
          test_vectors(&y_sp.vector, &y_gsl.vector, 1.0e-9,
                       "test_dgemv: compressed column format");

          gsl_spmatrix_free(Sc);
        }
    }

  gsl_matrix_free(A);
  gsl_vector_free(x);
  gsl_vector_free(y1);
  gsl_vector_free(y2);
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
