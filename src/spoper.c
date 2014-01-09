/* spoper.c
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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

#include "gsl_spmatrix.h"

int
gsl_spmatrix_scale(gsl_spmatrix *m, const double x)
{
  size_t i;

  for (i = 0; i < m->nz; ++i)
    m->data[i] *= x;

  return GSL_SUCCESS;
} /* gsl_spmatrix_scale() */

int
gsl_spmatrix_minmax(const gsl_spmatrix *m, double *min_out, double *max_out)
{
  double min, max;
  size_t n;

  if (m->nz == 0)
    {
      GSL_ERROR("matrix is empty", GSL_EINVAL);
    }

  min = m->data[0];
  max = m->data[0];

  for (n = 1; n < m->nz; ++n)
    {
      double x = m->data[n];

      if (x < min)
        min = x;

      if (x > max)
        max = x;
    }

  *min_out = min;
  *max_out = max;

  return GSL_SUCCESS;
} /* gsl_spmatrix_minmax() */

/*
gsl_spmatrix_d2sp()
  Convert a dense gsl_matrix to sparse (triplet) format

Inputs: S - (output) sparse matrix in triplet format
        A - (input) dense matrix to convert
*/

int
gsl_spmatrix_d2sp(gsl_spmatrix *S, const gsl_matrix *A)
{
  int s = GSL_SUCCESS;
  size_t i, j;

  gsl_spmatrix_reset(S);

  for (i = 0; i < A->size1; ++i)
    {
      for (j = 0; j < A->size2; ++j)
        {
          double x = gsl_matrix_get(A, i, j);

          if (x != 0.0)
            gsl_spmatrix_set(S, i, j, x);
        }
    }

  return s;
} /* gsl_spmatrix_d2sp() */

/*
gsl_spmatrix_sp2d()
  Convert a sparse matrix to dense format
*/

int
gsl_spmatrix_sp2d(gsl_matrix *A, const gsl_spmatrix *S)
{
  int s = GSL_SUCCESS;

  gsl_matrix_set_zero(A);

  if (GSLSP_ISTRIPLET(S))
    {
      size_t n;

      for (n = 0; n < S->nz; ++n)
        {
          size_t i = S->i[n];
          size_t j = S->j[n];
          double x = S->data[n];

          gsl_matrix_set(A, i, j, x);
        }
    }
  else
    {
      GSL_ERROR("non-triplet formats not yet supported", GSL_EINVAL);
    }

  return s;
} /* gsl_spmatrix_sp2d() */
