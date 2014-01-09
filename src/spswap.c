/* spswap.c
 * 
 * Copyright (C) 2014 Patrick Alken
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
gsl_spmatrix_transpose_memcpy(gsl_spmatrix *dest, const gsl_spmatrix *src)
{
  dest->size1 = src->size2;
  dest->size2 = src->size1;
  dest->nz = src->nz;
  dest->flags = src->flags;

  /* enlarge dest if needed */
  if (dest->nzmax < src->nzmax)
    {
      int s = gsl_spmatrix_realloc(src->nzmax, dest);
      if (s)
        return s;
    }

  if (GSLSP_ISTRIPLET(src))
    {
      size_t n;

      for (n = 0; n < src->nz; ++n)
        {
          dest->i[n] = src->p[n];
          dest->p[n] = src->i[n];
          dest->data[n] = src->data[n];
        }
    }
  else
    {
      GSL_ERROR("non-triplet formats not yet supported", GSL_EINVAL);
    }

  return GSL_SUCCESS;
} /* gsl_spmatrix_transpose_memcpy() */
