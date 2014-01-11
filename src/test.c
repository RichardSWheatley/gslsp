/* test.c
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
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>

#include "gsl_spmatrix.h"
#include "test.h"

/*
create_random_sparse()
  Create a random sparse matrix with approximately
M*N*density non-zero entries

Inputs: M       - number of rows
        N       - number of columns
        density - sparse density \in (0,1)
                  0 = no non-zero entries
                  1 = all m*n entries are filled
        r       - random number generator

Return: pointer to sparse matrix in triplet format (must be freed by caller)
*/

static gsl_spmatrix *
create_random_sparse(const size_t M, const size_t N, const double density,
                     const gsl_rng *r)
{
  gsl_spmatrix *m = gsl_spmatrix_alloc(M, N);
  size_t nnzwanted = (size_t) round(M * N * GSL_MIN(density, 1.0));
  size_t n = 0;

  while (n <= nnzwanted)
    {
      /* generate a random row and column */
      size_t i = gsl_rng_uniform(r) * M;
      size_t j = gsl_rng_uniform(r) * N;
      double x;

      assert(i < M);
      assert(j < N);

      /* check if this row/column is already filled */
      if (gsl_spmatrix_get(m, i, j) != 0.0)
        continue;

      /* generate random m_{ij} and add it */
      x = gsl_rng_uniform(r);
      gsl_spmatrix_set(m, i, j, x);
      ++n;
    }

  return m;
} /* create_random_sparse() */

void
test_random_sparse_matrix(gsl_matrix *m, gsl_rng *r, double lower,
                          double upper)
{
  size_t M = m->size1;
  size_t N = m->size2;
  size_t i, j;
  int numel; /* number of non-zero elements in a row */
  double x;

  if (N == 1 && M == 1)
    {
      x = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_matrix_set(m, 0, 0, x);
      return;
    }

  gsl_matrix_set_zero(m);

  for (i = 0; i < M; ++i)
    {
      /* pick a random number between 1 and N/2 - this is how many
       * nonzero elements are in this row
       */
      numel = (int) (gsl_rng_uniform(r) * (N / 2) + 1);
      for (j = 0; j < (size_t) numel; ++j)
        {
          int k = (int) (gsl_rng_uniform(r) * (N - 1)); /* pick random column in [0,N-1] */
          x = gsl_rng_uniform(r) * (upper - lower) + lower;
          gsl_matrix_set(m, i, k, x);
        }

      /* always set the diagonal element if possible */
      if (i < N)
        {
          x = gsl_rng_uniform(r) * (upper - lower) + lower;
          gsl_matrix_set(m, i, i, x);
        }
    }
} /* test_random_sparse_matrix() */

void
test_random_vector(gsl_vector *v, gsl_rng *r, double lower,
                   double upper)
{
  size_t N = v->size;
  size_t i;
  double x;

  for (i = 0; i < N; ++i)
    {
      x = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_vector_set(v, i, x);
    }
} /* test_random_vector() */

int
test_vectors(gsl_vector *observed, gsl_vector *expected, const double tol,
             const char *str)
{
  int s = 0;
  size_t N = observed->size;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      double x_obs = gsl_vector_get(observed, i);
      double x_exp = gsl_vector_get(expected, i);

      gsl_test_rel(x_obs, x_exp, tol, "N=%zu i=%zu %s", N, i, str);
    }

  return s;
} /* test_vectors() */

static void
test_getset(const size_t M, const size_t N, const gsl_rng *r)
{
  int status;
  size_t i, j;
  size_t k = 0;
  gsl_spmatrix *m = gsl_spmatrix_alloc(M, N);

  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double x = (double) ++k;
          double y;

          gsl_spmatrix_set(m, i, j, x);
          y = gsl_spmatrix_get(m, i, j);

          if (x != y)
            status = 1;
        }
    }

  gsl_test(status, "test_getset: M=%zu N=%zu _get != _set", M, N);

  gsl_spmatrix_free(m);
} /* test_getset() */

static void
test_memcpy(const size_t M, const size_t N, const gsl_rng *r)
{
  int status;
  gsl_spmatrix *a = create_random_sparse(M, N, 0.2, r);
  gsl_spmatrix *b, *c, *d;
  
  b = gsl_spmatrix_memcpy(a);

  status = gsl_spmatrix_equal(a, b) != 1;
  gsl_test(status, "test_memcpy: M=%zu N=%zu triplet format", M, N);

  c = gsl_spmatrix_compcol(a);
  d = gsl_spmatrix_memcpy(c);

  status = gsl_spmatrix_equal(c, d) != 1;
  gsl_test(status, "test_memcpy: M=%zu N=%zu compressed column format", M, N);

  gsl_spmatrix_free(a);
  gsl_spmatrix_free(b);
  gsl_spmatrix_free(c);
  gsl_spmatrix_free(d);
} /* test_memcpy() */

int
main()
{
  const size_t N_max = 20;
  size_t M, N;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  for (M = 1; M <= N_max; ++M)
    {
      for (N = 1; N <= N_max; ++N)
        {
          test_memcpy(M, N, r);
          test_getset(M, N, r);
        }
    }

#if 0
  fprintf(stderr, "test: getset...");
  s += test_getset();
  fprintf(stderr, "done (s = %d)\n", s);

  fprintf(stderr, "test: dgemv...");
  s += test_dgemv();
  fprintf(stderr, "done (s = %d)\n", s);
#endif

  gsl_rng_free(r);

  exit (gsl_test_summary());
} /* main() */
