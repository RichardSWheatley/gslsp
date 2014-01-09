/* spcompress.c
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

#include "gsl_spmatrix.h"

static int cumulative_sum(size_t *c, size_t *p, size_t n);

/*
gsl_spmatrix_compcol()
  Create a sparse matrix in compressed column format

Inputs: T - sparse matrix in triplet format

Return: pointer to new matrix (should be freed when finished with it)
*/

gsl_spmatrix *
gsl_spmatrix_compcol(const gsl_spmatrix *T)
{
  const size_t *Tj = T->p;
  gsl_spmatrix *m;
  size_t *c; /* counts of non-zeros in each column */
  size_t nz = T->nz;
  size_t n;

  m = gsl_spmatrix_alloc_nzmax(T->size1, T->size2, nz,
                               GSL_SPMATRIX_COMPCOL);
  if (!m)
    return NULL;

  /* allocate column pointer array */
  c = calloc(1, (T->size2 + 1) * sizeof(size_t));
  if (!c)
    {
      GSL_ERROR_VAL("failed to allocate column pointers",
                    GSL_ENOMEM, 0);
    }

  /*
   * compute the number of elements in each column:
   * c[j] = # non-zero elements in column j
   */
  for (n = 0; n < nz; ++n)
    c[Tj[n]]++;

  /*
   * now convert to column pointers:
   *
   * p[j] = Sum_{i=0...j-1} c[i]
   */
  cumulative_sum(c, m->p, m->size2);

  for (n = 0; n < nz; ++n)
    {
      size_t k = c[Tj[n]]++;

      m->i[k] = T->i[n];
      m->data[k] = T->data[n];
    }

  free(c);

  m->nz = T->nz;
  m->flags = GSL_SPMATRIX_COMPCOL;

  return m;
} /* gsl_spmatrix_compcol() */

/*
cumulative_sum
  Compute an in-place cumulative sum. Given an array c, compute:

p[i] = Sum_{j=0...i-1} c[j]

or as a recurrence relation:

p[0] = 0
p[i] = p[i-1] + c[i-1]

Inputs: c - input array c
        p - (output) output array of cumulative sums
        n - number of elements in c array on input

Notes:

1) output array p is of size (n + 1), with p[n] = Sum_j c[j]

2) on output, p[0:n-1] is copied to c[0:n-1]
*/

static int
cumulative_sum(size_t *c, size_t *p, size_t n)
{
  int s = GSL_SUCCESS;
  size_t sum = 0;
  size_t k;

  for (k = 0; k < n; ++k)
    {
      p[k] = sum;
      sum += c[k];

      c[k] = p[k];
    }

  p[n] = sum;

  return s;
} /* cumulative_sum() */
