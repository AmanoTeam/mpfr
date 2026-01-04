/* Test file for (physicist's) hermite polynomials.

Copyright 2025 Free Software Foundation, Inc.
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

#define RANDOM_TESTS_BATCH 50

static void
test_special_cases (void)
{
  mpfr_t x, expected, res;
  int ret, rnd, degree_2 = 2, degree_4 = 4, odd_degree = 3;

  mpfr_init2 (res, 100);
  mpfr_init2 (expected, 200);

  /* init x = 0 */
  mpfr_init2 (x, 10);
  MPFR_SET_ZERO (x);

  RND_LOOP (rnd)
    {
      mpfr_rnd_t r = (mpfr_rnd_t) rnd;
      if (mpfr_hermite (res, odd_degree, x, r) != 0
          || !MPFR_IS_ZERO (res))
        {
          printf ("H_%d(0) should be 0; got ", odd_degree);
          mpfr_out_str (stdout, 10, 0, res, r);
          printf ("\n");
          exit (1);
        }

      /* H_2(0) = -2 */
      mpfr_set_d (expected, -2.0, r);
      ret = mpfr_hermite (res, degree_2, x, r);
      if (ret != 0 || !mpfr_equal_p (res, expected))
        {
          printf ("H_2(0) should be: ");
          mpfr_dump (expected);
          printf ("got: ");
          mpfr_dump (res);
          printf ("with return value: %d\n", ret);
          exit (1);
        }

      /* H_4(0) = 12 */
      mpfr_set_d (expected, 12.0, r);
      ret = mpfr_hermite (res, degree_4, x, r);
      if (ret != 0 || !mpfr_equal_p (res, expected))
        {
          printf ("H_4(0) should be: ");
          mpfr_dump (expected);
          printf ("got: ");
          mpfr_dump (res);
          printf ("with return value: %d\n", ret);
          exit (1);
        }

      /* H_{1e6}(0) exceed 2^13, therefore should return NaN */
      ret = mpfr_hermite (res, (unsigned) 1e6, x, r);
      if (ret != 0 || !mpfr_nan_p (res))
        {
          printf ("H_{1e6}(0) return NaN. Got: ");
          mpfr_dump (res);
          printf ("with return value: %d\n", ret);
          exit (1);
        }
    }

  mpfr_clear (x);
  mpfr_clear (res);
  mpfr_clear (expected);
  mpfr_free_cache ();
}

/* for n == 0:
     - H_0(Inf)  -> 1.0
     - H_0(-Inf) -> 1.0
     - H_0(NaN)  -> 1.0
   for n >= 1:
     - H_n(Inf)  -> NaN
     - H_n(-Inf) -> NaN
     - H_n(NaN)  -> NaN */
static void
test_singular_input (void)
{
  mpfr_t res, x;
  int i, ret;
  unsigned n, degrees[] = { 0, 1, 2, 5, 12, 21, 30, 50 };

  mpfr_init2 (res, 200);
  mpfr_init2 (x, 200);

  for (i = 0; i < sizeof (degrees) / sizeof (unsigned); i++)
    {
      n = degrees[i];

      mpfr_set_nan (x);
      ret = mpfr_hermite (res, n, x, MPFR_RNDN);
      if (ret != 0 || !mpfr_nan_p (res))
        {
          printf ("For x = NAN, H_%u should be NAN\ngot: ", n);
          mpfr_dump (res);
          printf ("With return value: %d\n", ret);
          exit (1);
        }
      mpfr_set_inf(x, 1);
      ret = mpfr_hermite (res, n, x, MPFR_RNDN);
      if (ret != 0 || !mpfr_nan_p (res))
        {
          printf ("For x = +Inf, H_%u should be NAN\ngot: ", n);
          mpfr_dump (res);
          printf ("With return value: %d\n", ret);
          exit (1);
        }
      mpfr_set_inf(x, -1);
      ret = mpfr_hermite (res, n, x, MPFR_RNDN);
      if (ret != 0 || !mpfr_nan_p (res))
        {
          printf ("For x = -Inf, H_%u should be NAN\ngot: ", n);
          mpfr_dump (res);
          printf ("With return value: %d\n", ret);
          exit (1);
        }
    }

  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_free_cache ();
}

/* hermite(0,x) -> 1.0 */
static void
test_first_iteration (void)
{
  mpfr_t one, res, x;
  int ret;

  mpfr_init2 (res, 200);
  mpfr_init2 (x, 200);

  mpfr_set_d (x, 0.94, MPFR_RNDD);

  /* init one constant */
  mpfr_init2 (one, 200);
  mpfr_set_ui (one, 1, MPFR_RNDD);

  /* The first iteration should always be 1 */
  ret = mpfr_hermite (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_equal_p (res, one))
    {
      printf ("The first Hermite polynomial H_0 should be exactly 1.\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }

  mpfr_clear (one);
  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_free_cache ();
}

/* hermite(1,x) -> 2x */
static void
test_second_iteration (void)
{
  mpfr_t x, expected, res;
  int ret;

  mpfr_init2 (x, 100);
  mpfr_init2 (res, 100);
  mpfr_init2 (expected, 200);

  mpfr_set_d (x, 1.0/3.0, MPFR_RNDD);
  mpfr_mul_ui (expected, x, 2, MPFR_RNDN);

  /* The second iteration should always return 2x. Since prec(res) = prec(x),
     ret should be 0 */
  ret = mpfr_hermite (res, 1, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_equal_p (res, expected))
    {
      printf ("H_1 should be 2x\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }

  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_clear (expected);
  mpfr_free_cache ();
}

static void
test_double_precision (void)
{
  mpfr_t x, expected, res;

  mpfr_init2 (x, IEEE754_DOUBLE_PREC);
  mpfr_init2 (res, IEEE754_DOUBLE_PREC);
  mpfr_init2 (expected, IEEE754_DOUBLE_PREC);

  /* H_3(3.49376) = 2.9924358881463501e2 with MPFR_RNDN */
  mpfr_set_d (x, 3.49376, MPFR_RNDD);
  mpfr_set_str (expected, "299.24358881463500799999999999999999999961", 10, MPFR_RNDN);
  mpfr_hermite (res, 3, x, MPFR_RNDN);
  if (mpfr_cmp (res, expected))
    {
      printf ("test_double_precision failed;\n");
      DUMP_NUMBERS (expected, res);
      exit (1);
    }

  /* H_6(-2.2364) = -5.1892977013945506e2 with MPFR_RNDN */
  mpfr_set_d (x, -2.2364, MPFR_RNDD);
  mpfr_set_str (expected, "-518.92977013945504945520448445364119704459", 10, MPFR_RNDN);
  mpfr_hermite (res, 6, x, MPFR_RNDN);
  if (mpfr_cmp (res, expected))
    {
      printf ("test_double_precision failed;\n");
      DUMP_NUMBERS (expected, res);
      exit (1);
    }

  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_clear (expected);
  mpfr_free_cache ();
}

int
main (void)
{
  tests_start_mpfr ();

  /* tests for H_n(0) */
  test_special_cases ();

  /* for singular input (+/-Inf or NaN), mpfr_hermit should always set res
     to NaN and return 0 */
  test_singular_input ();

  /* the first two iterations are tested separately because they are the base
     cases of the recursion algorithm used to calculate mpfr_hermite */
  test_first_iteration ();
  test_second_iteration ();

  test_double_precision ();

  random_poly_suite (RANDOM_TESTS_BATCH, RANDOM_TESTS_BATCH,
                     IEEE754_DOUBLE_PREC, mpfr_legendre);

  tests_end_mpfr ();
  return 0;
}
