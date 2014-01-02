/* test.h
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

#ifndef __TEST_H__
#define __TEST_H__

void test_random_sparse_matrix(gsl_matrix *m, gsl_rng *r, double lower,
                               double upper);
void test_random_vector(gsl_vector *v, gsl_rng *r, double lower,
                        double upper);
int test_vectors(gsl_vector *observed, gsl_vector *expected,
                 const double tol, const char *str);

int test_getset(void);
int test_dgemv(void);

#endif /* __TEST_H__ */
