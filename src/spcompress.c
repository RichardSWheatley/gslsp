/* spcompress.c
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include "gsl_spmatrix.h"

static int cumulative_sum(const size_t n, size_t *c);

/*
gsl_spmatrix_compcol()
  Create a sparse matrix in compressed column format

Inputs: T - sparse matrix in triplet format

Return: pointer to new matrix (should be freed when finished with it)
*/

gsl_spmatrix *
gsl_spmatrix_compcol(const gsl_spmatrix *T)
{
  const size_t *Tj; /* column indices of triplet matrix */
  size_t *Cp;       /* column pointers of compressed column matrix */
  gsl_spmatrix *m;
  size_t *c; /* counts of non-zeros in each column */
  size_t nz = T->nz;
  size_t n;

  m = gsl_spmatrix_alloc_nzmax(T->size1, T->size2, nz,
                               GSL_SPMATRIX_COMPCOL);
  if (!m)
    return NULL;

  /* allocate column pointer array */
  c = calloc(1, T->size2 * sizeof(size_t));
  if (!c)
    {
      GSL_ERROR_VAL("failed to allocate workspace",
                    GSL_ENOMEM, 0);
    }

  Tj = T->p;
  Cp = m->p;

  /* initialize column pointers to 0 */
  for (n = 0; n < m->size2 + 1; ++n)
    Cp[n] = 0;

  /*
   * compute the number of elements in each column:
   * Cp[j] = # non-zero elements in column j
   */
  for (n = 0; n < nz; ++n)
    Cp[Tj[n]]++;

  /* compute column pointers: p[j] = p[j-1] + nnz[j-1] */
  cumulative_sum(m->size2, Cp);

  /* make a copy of the column pointers */
  for (n = 0; n < m->size2; ++n)
    c[n] = Cp[n];

  /* transfer data from triplet format to compressed column */
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
cumulative_sum()

Compute the cumulative sum:

p[j] = Sum_{k=0...j-1} c[k]

0 <= j < n + 1

Alternatively,
p[0] = 0
p[j] = p[j - 1] + c[j - 1]

Inputs: n - length of input array
        c - (input/output) array of size n + 1
            on input, contains the n values c[k]
            on output, contains the n + 1 values p[j]

Return: success or error
*/

static int
cumulative_sum(const size_t n, size_t *c)
{
  size_t sum = 0;
  size_t k;

  for (k = 0; k < n; ++k)
    {
      size_t ck = c[k];
      c[k] = sum;
      sum += ck;
    }

  c[n] = sum;

  return GSL_SUCCESS;
} /* cumulative_sum() */
