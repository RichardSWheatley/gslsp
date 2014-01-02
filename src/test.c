/* test.c
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>

#include "gsl_spmatrix.h"
#include "test.h"

void
test_random_sparse_matrix(gsl_matrix *m, gsl_rng *r, double lower,
                          double upper)
{
  size_t N = m->size1;
  size_t M = m->size2;
  size_t i, j;
  int numel; /* number of non-zero elements in a row */
  double x;

  if (N == 1)
    {
      x = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_matrix_set(m, 0, 0, x);
      return;
    }

  gsl_matrix_set_zero(m);

  for (i = 0; i < N; ++i)
    {
      /* pick a random number between 1 and M/2 - this is how many
       * nonzero elements are in this row
       */
      numel = (int) (gsl_rng_uniform(r) * (M / 2 - 1) + 1);
      for (j = 0; j < (size_t) numel; ++j)
        {
          int k = (int) (gsl_rng_uniform(r) * (M - 2));
          x = gsl_rng_uniform(r) * (upper - lower) + lower;
          gsl_matrix_set(m, i, k, x);
        }

      /* always set the diagonal element */
      x = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_matrix_set(m, i, i, x);
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

      gsl_test_rel(x_obs, x_exp, tol, str);
    }

  return s;
} /* test_vectors() */

int
main()
{
  int s = 0;

  fprintf(stderr, "test: getset...");
  s += test_getset();
  fprintf(stderr, "done (s = %d)\n", s);

  fprintf(stderr, "test: dgemv...");
  s += test_dgemv();
  fprintf(stderr, "done (s = %d)\n", s);

  exit (gsl_test_summary());
} /* main() */
