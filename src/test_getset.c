/* test_getset.c
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

#include "gsl_spmatrix.h"
#include "test.h"

int
test_getset(void)
{
  int s = 0;
  size_t N_max = 50;
  size_t N;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_spmatrix *S = gsl_spmatrix_alloc(1, 1);

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix *A = gsl_matrix_calloc(N, N);
      size_t i, j;

      /* create random matrix */
      test_random_sparse_matrix(A, r, -10.0, 10.0);

      gsl_spmatrix_reset(S);

      for (i = 0; i < N; ++i)
        {
          for (j = 0; j < N; ++j)
            {
              double x = gsl_matrix_get(A, i, j);
              double y;

              gsl_spmatrix_set(S, i, j, x);
              y = gsl_spmatrix_get(S, i, j);

              if (x != 0.0)
                {
                  gsl_test_rel(y, x, 1.0e-10, "N = %zu, i = %zu, j = %zu",
                               N, i, j);
                }
              else
                {
                  gsl_test_abs(y, 0.0, 1.0e-10, "N = %zu, i = %zu, j = %zu",
                               N, i, j);
                }
            }
        }

      gsl_matrix_free(A);
    }

  gsl_rng_free(r);
  gsl_spmatrix_free(S);

  return s;
} /* test_getset() */
