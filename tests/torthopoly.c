/* Random test suite for  orthogonal polynomials.

Copyright 2025-2026 Free Software Foundation, Inc.
Contributed by Matteo Nicoli.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

#include "mpfr-test.h"

#ifndef X_LOWER_BOUND
#define X_LOWER_BOUND -1.0
#endif

#ifndef X_HIGHER_BOUND
#define X_HIGHER_BOUND 1.0
#endif

static void
random_array (int *array, int size, int inf, int sup)
{
  int i, tmp;
  for (i = 0; i < size; i++)
    {
      if (inf > sup)
        {
          tmp = inf;
          inf = sup;
          sup = tmp;
        }
      array[i] = inf + (int) (randulong () % (unsigned long) (sup - inf + 1));
    }
}

static void
test_poly_random (int n, mpfr_prec_t p, unsigned long K)
{
  mpfr_t x, y, z, t;
  unsigned long k;
  int rnd;
  mpfr_init2 (x, p);
  mpfr_init2 (y, p);
  mpfr_init2 (z, p + 20);
  mpfr_init2 (t, p);
  for (k = 0; k < K; k++)
    {
      mpfr_urandomb (x, RANDS); /* x is in [0,1] */
      mpfr_mul_d (x, x, X_HIGHER_BOUND - X_LOWER_BOUND, MPFR_RNDN);
      mpfr_add_d (x, x, X_LOWER_BOUND, MPFR_RNDN);
      /* now x is in [X_LOWER_BOUND, X_HIGHER_BOUND] */
      RND_LOOP_NO_RNDF (rnd)
        {
          MPFR_ORTHOGONAL_POLY_FN (y, n, x, (mpfr_rnd_t) rnd);
          MPFR_ORTHOGONAL_POLY_FN (z, n, x, MPFR_RNDN);
          if (mpfr_can_round (z, p + 20, MPFR_RNDN, (mpfr_rnd_t) rnd, p))
            {
              mpfr_set (t, z, (mpfr_rnd_t) rnd);
              if (mpfr_cmp (y, t))
                {
                  printf ("Error in random test for n=%d x=", n);
                  mpfr_out_str (stdout, 16, 0, x, MPFR_RNDN);
                  printf (" rnd=%s\n", mpfr_print_rnd_mode ((mpfr_rnd_t) rnd));
                  DUMP_NUMBERS (t, y);
                  exit (1);
                }
            }
        }
    }
  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);
}

static void
random_poly_suite (int num_degrees, int ntests, mpfr_prec_t p)
{
  /* we set the minimum degree to 2 to skip the two base cases P0 and P1,
     and the maximum degree to 128 to limit the range of degrees tested
     to the same limit of the C++ standard */
  int i, min_degree = 2, max_degree = 128;
  int *test_degrees;

  test_degrees = (int *) malloc (num_degrees * sizeof(int));
  if (!test_degrees)
    {
      printf ("Could not allocate memory for random tests\n");
      exit (1);
    }

  random_array (test_degrees, num_degrees, min_degree, max_degree);

  for (i = 0; i < num_degrees; i++)
    {
      test_poly_random (test_degrees[i], p, ntests);
    }

  free (test_degrees);
}
