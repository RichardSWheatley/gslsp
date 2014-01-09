/* spgetset.c
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

double
gsl_spmatrix_get(const gsl_spmatrix *m, const size_t i, const size_t j)
{
  if (!(m->flags & GSL_SPMATRIX_TRIPLET))
    {
      GSL_ERROR("matrix not in triplet representation", GSL_EINVAL);
    }
  else
    {
      size_t n;
      const size_t *Ti = m->i;
      const size_t *Tj = m->p;

      for (n = 0; n < m->nz; ++n)
        {
          if (Ti[n] == i && Tj[n] == j)
            return m->data[n];
        }

      /* element not found; return 0 */
      return 0.0;
    }
} /* gsl_spmatrix_get() */

/*
gsl_spmatrix_set()
  Add an element to a matrix in triplet form

Inputs: m - spmatrix
        i - row index
        j - column index
        x - matrix value
*/

int
gsl_spmatrix_set(gsl_spmatrix *m, const size_t i, const size_t j,
                 const double x)
{
  if (!(m->flags & GSL_SPMATRIX_TRIPLET))
    {
      GSL_ERROR("matrix not in triplet representation", GSL_EINVAL);
    }
  else if (x == 0.0)
    return GSL_SUCCESS;
  else
    {
      int s = GSL_SUCCESS;

      if (m->nz >= m->nzmax)
        {
          s = gsl_spmatrix_realloc(2 * m->nzmax, m);
          if (s)
            return s;
        }

      /* store the triplet (i, j, x) */
      m->i[m->nz] = i;
      m->p[m->nz] = j;
      m->data[m->nz] = x;

      /* increase matrix dimensions if needed */
      m->size1 = GSL_MAX(m->size1, i + 1);
      m->size2 = GSL_MAX(m->size2, j + 1);

      ++(m->nz);

      return s;
    }
} /* gsl_spmatrix_set() */
