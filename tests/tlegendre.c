/* Test file for legendre polynomials.

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

#define ARBITRARILY_LOW_PREC 10
#define IEEE754_SINGLE_PREC  24
#define IEEE754_DOUBLE_PREC  53
#define MPFR_PREC_100        100
#define MPFR_PREC_200        200

static const unsigned degrees[] =
{
  2,       /* The first even degree after the base cases */
  3,       /* The first odd degree after the base cases */
  10,
  50,
  128,     /* The maximum degree officially supported by the C++ standard */
  1024,    /* 2^10 */
  8192,    /* 2^13 */
  1048576, /* 2^20 */
};
static const char *expected_vals[] =
{
  "-0.00100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
  "-0.01110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
  "-0.00110000001011111100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
  "-0.11111110011011111010011011101111010001100110010000011011100111001011011101000001101100101111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e-5",
  "-0.10100000000001110010100100000001010100100100100110100000000110011010001001111001010101100101000100110011011001000100100000001011011110011111100010011110111000001000111000111001110011111111101000111101e-5",
  "-0.10011011001011010000011110100001001010000100011101101001111011101111100001001000000001001001111000111010100010110101101100101110100011000100010001101010010001111010101001110011010101011000110011001001e-5",
  "-0.10100000101010101110101000101101101011000101000000110110011000110011100011010000111011000101100110010101010100011110011001001111001111110111101101010000100110001110101100110000010110100001001110000100e-8",
  "-0.10011011001100111110110010111111011011111101111001101010111101000000001111110101001101100011100101011001000110111010100000110000101101111101100001001111100111001111111011111101111001000000011000100010e-10",
};

static const int n_degrees_test = sizeof(degrees)/sizeof(unsigned);

static void
test_domain (void)
{
  mpfr_t res, upper, lower, inner, outer;

  mpfr_init2 (res, 200);

  /* init the constants to (respectively):
       -  the upper bound of the domain
       - the lower bound of the domain
       - a number within the dmain
       - a nomber outside of the domain */
  mpfr_init2 (upper, 10);
  mpfr_set_ui (upper, 1, MPFR_RNDD);
  mpfr_init2 (lower, 200);
  mpfr_set_d (lower, -1.0, MPFR_RNDD);
  mpfr_init2 (inner, 200);
  mpfr_set_d (inner, 1.0/5.0, MPFR_RNDD);
  mpfr_init2 (outer, 200);
  /* 1.00000000001 */
  mpfr_set_d (outer, 1.0 + 1e-11, MPFR_RNDD);

  int i;
  for (i = 0; i < 10; i++)
    {
      mpfr_legendre (res, i, upper, MPFR_RNDN);
      if (MPFR_IS_NAN (res))
        {
          fprintf (stderr, "Upper bound input value ");
          mpfr_out_str (stderr, 10, 0, upper, MPFR_RNDD);
          fprintf (stderr, " should *not* lead to a NAN result (degree: %d)\n", i);
          exit (1);
        }
    }

  for (i = 0; i < 10; i++)
    {
      mpfr_legendre (res, i, lower, MPFR_RNDN);
      if (MPFR_IS_NAN (res))
        {
          fprintf (stderr, "Lower bound input value ");
          mpfr_out_str (stderr, 10, 0, lower, MPFR_RNDD);
          fprintf (stderr, " should *not* lead to a NAN result (degree: %d)\n", i);
          exit (1);
        }
    }

  for (i = 0; i < 10; i++)
    {
      mpfr_legendre (res, i, inner, MPFR_RNDN);
      if (MPFR_IS_NAN (res))
        {
          fprintf (stderr, "input number ");
          mpfr_out_str (stderr, 10, 0, inner, MPFR_RNDD);
          fprintf (stderr, " should *not* lead to a NAN result (degree: %d)\n", i);
          exit (1);
        }
    }

  for (i = 0; i < 10; i++)
    {
      mpfr_legendre (res, i, outer, MPFR_RNDN);
      if (!MPFR_IS_NAN (res))
        {
          fprintf (stderr, "input number ouside of the domain ");
          mpfr_out_str (stderr, 10, 0, outer, MPFR_RNDD);
          fprintf (stderr, " should lead to a NAN result (degree: %d)\n", i);
          exit (1);
        }
    }

  mpfr_clear (res);
  mpfr_clear (upper);
  mpfr_clear (lower);
  mpfr_clear (inner);
  mpfr_clear (outer);
  mpfr_free_cache ();
}

static void
test_domain_bounds (void)
{
  mpfr_t one, minus_one, res;
  int even_degree = 2, odd_degree = 3;

  mpfr_init2 (res, 200);

  /* init 1 constant */
  mpfr_init2 (one, 200);
  mpfr_set_ui (one, 1, MPFR_RNDD);

  /* init -1 constant */
  mpfr_init2 (minus_one, 200);
  mpfr_set_si (minus_one, -1, MPFR_RNDD);

  if (mpfr_legendre (res, even_degree, one, MPFR_RNDD) != 0
      || !mpfr_equal_p (one, res))
    {
      fprintf (stderr, "P%d(1) should be 1; got ", even_degree);
      mpfr_out_str (stderr, 10, 0, res, MPFR_RNDD);
      fprintf (stderr, "\n");
      exit (1);
    }

  if (mpfr_legendre (res, odd_degree, one, MPFR_RNDD) != 0
      || !mpfr_equal_p (one, res))
    {
      fprintf (stderr, "P%d(1) should be 1; got ", odd_degree);
      mpfr_out_str (stderr, 10, 0, res, MPFR_RNDD);
      fprintf (stderr, "\n");
      exit (1);
    }

  if (mpfr_legendre (res, even_degree, minus_one, MPFR_RNDD) != 0
      || !mpfr_equal_p (one, res))
    {
      fprintf (stderr, "P%d(-1) should be 1; got ", even_degree);
      mpfr_out_str (stderr, 10, 0, res, MPFR_RNDD);
      fprintf (stderr, "\n");
      exit (1);
    }

  if (mpfr_legendre (res, odd_degree, minus_one, MPFR_RNDD) != 0
      || !mpfr_equal_p (minus_one, res))
    {
      fprintf (stderr, "P%d(-1) should be -1; got ", odd_degree);
      mpfr_out_str (stderr, 10, 0, res, MPFR_RNDD);
      fprintf (stderr, "\n");
      exit (1);
    }

  mpfr_clear (one);
  mpfr_clear (minus_one);
  mpfr_clear (res);
  mpfr_free_cache ();
}

/* legendre(0,x) -> 1.0 */
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

  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  /* The first iteration should always be 1 */
  if (ret != 0 || !mpfr_equal_p (res, one))
    {
      fprintf (stderr, "The first legendre polynomial P_0 should be exactly 1.\ngot: ");
      mpfr_dump (res);
      fprintf (stderr, "With return value: %d\n", ret);
      exit (1);
    }

  /* The first iteration should ignore x, but still check
     that x is within the canonical domain [-1,1]*/
  mpfr_set_zero (x, 1);
  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_equal_p (res, one))
    {
      fprintf (stderr, "P_0 should be 1, if x = +0.0\n got: ");
      mpfr_dump (res);
      fprintf (stderr, "With return value: %d\n", ret);
      exit (1);
    }
  mpfr_set_zero (x, -1);
  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_equal_p (res, one))
    {
      fprintf (stderr, "P_0 should be 1, if x = -0.0\ngot: ");
      mpfr_dump (res);
      fprintf (stderr, "With return value: %d\n", ret);
      exit (1);
    }

  /* For x equal to NAN, +Inf or -Inf, the result of the first iteration
     should be NAN */
  mpfr_set_nan (x);
  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_nan_p (res))
    {
      fprintf (stderr, "For x = NAN, P_0 should be NAN\ngot: ");
      mpfr_dump (res);
      fprintf (stderr, "With return value: %d\n", ret);
      exit (1);
    }
  mpfr_set_inf(x, 1);
  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_nan_p (res))
    {
      fprintf (stderr, "For x = +Inf, P_0 should be NAN\ngot: ");
      mpfr_dump (res);
      fprintf (stderr, "With return value: %d\n", ret);
      exit (1);
    }
  mpfr_set_inf(x, -1);
  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_nan_p (res))
    {
      fprintf (stderr, "For x = -Inf, P_0 should be NAN\ngot: ");
      mpfr_dump (res);
      fprintf (stderr, "With return value: %d\n", ret);
      exit (1);
    }

  mpfr_clear (one);
  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_free_cache ();
}

/* legendre(1,x) -> x */
static void
test_second_iteration (void)
{
  mpfr_t res, x;
  int ret;

  mpfr_init2 (res, 200);
  mpfr_init2 (x, 200);

  mpfr_set_d (x, 1.0/3.0, MPFR_RNDD);

  /* The second iteration should return x */
  ret = mpfr_legendre (res, 1, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_equal_p (res, x))
    {
      fprintf (stderr, "P_1 should be x itself\ngot: ");
      mpfr_dump (res);
      fprintf (stderr, "With return value: %d\n", ret);
      exit (1);
    }

  /* For x equal to NAN, +Inf or -Inf, the result of the second iteration
     should be NAN as well */
  mpfr_set_nan (x);
  ret = mpfr_legendre (res, 1, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_nan_p (res))
    {
      fprintf (stderr, "For x = NAN, P_1 should be NAN\n got: ");
      mpfr_dump (res);
      fprintf (stderr, "With return value: %d\n", ret);
      exit (1);
    }
  mpfr_set_inf(x, 1);
  ret = mpfr_legendre (res, 1, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_nan_p (res))
    {
      fprintf (stderr, "For x = +Inf, P_1 should be NAN\ngot: ");
      mpfr_dump (res);
      fprintf (stderr, "With return value: %d\n", ret);
      exit (1);
    }
  mpfr_set_inf(x, -1);
  ret = mpfr_legendre (res, 1, x, MPFR_RNDN);
  if (ret !=0 || !mpfr_nan_p (res))
    {
      fprintf (stderr, "For x = -Inf, P_1 should be NAN\ngot: ");
      mpfr_dump (res);
      fprintf (stderr, "With return value: %d\n", ret);
      exit (1);
    }

  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_free_cache ();
}

static void
test_sample_with_precision (unsigned x_prec, unsigned res_prec)
{
  int ret;
  mpfr_t res, x, expected;
  const char *x_val = "0.5";

  mpfr_init2 (x, x_prec);
  mpfr_init2 (expected, res_prec);
  mpfr_init2 (res, res_prec);

  mpfr_set_str (x, x_val, 10, MPFR_RNDN);

  unsigned i = 0;
  for (; i < n_degrees_test; i++)
    {
      mpfr_set_str (expected, expected_vals[i], 2, MPFR_RNDN);
      ret = mpfr_legendre (res, degrees[i], x, MPFR_RNDN);

      if (ret != 0)
        {
          fprintf (stderr, "P%d(%s) should be exact. Got `%d` as ternary value\n",
                   degrees[i], x_val, ret);
          exit (1);
        }

      if (!mpfr_eq (res, expected,  mpfr_min_prec (expected)))
        {
          fprintf (stderr, "Wrong value for P%d(%s) [prec(x)=%u prec(res)=%u]\n",
                  degrees[i], x_val, x_prec, res_prec);
          fprintf (stderr, "expected ");
          mpfr_dump (expected);
          fprintf (stderr, "got      ");
          mpfr_dump (res);
          exit (1);
        }
    }

  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_clear (expected);
  mpfr_free_cache ();
}

static void
test_round (unsigned prec)
{
  const mpfr_rnd_t rounding_modes[] =
    {
      MPFR_RNDN,
      MPFR_RNDZ,
      MPFR_RNDU,
      MPFR_RNDD,
      MPFR_RNDA,
      MPFR_RNDF,
    };
  mpfr_t res, x, expected;
  const char *x_val = "0.5";
  int ret;

  mpfr_init2 (x, prec);
  mpfr_init2 (expected, prec);
  mpfr_init2 (res, prec);

  mpfr_set_str (x, x_val, 10, MPFR_RNDN);

  unsigned i = 0, j = 0;
  for (; j < sizeof(rounding_modes) / sizeof(int); j++)
    {
      for (; i < n_degrees_test; i++)
        {
          mpfr_set_str (expected, expected_vals[i], 2, MPFR_RNDN);
          ret = mpfr_legendre (res, degrees[i], x, rounding_modes[j]);

          if (ret != 0)
            {
              fprintf (stderr, "Wrong rounding [%s] for P%d(%s) [prec(x)=prec(res)=%u]\n",
                       mpfr_print_rnd_mode (rounding_modes[j]), degrees[i], x_val, prec);
              exit (1);
            }
        }
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

  /* The canonical domain of Legendre polynomials is [-1,1]. mpfr_legendre
     should return false for any x outside of the canonical domain */
  test_domain ();

  /* Upper and lower bounds are computed without using Bonnet's recursion */
  test_domain_bounds ();

  /* The first two iterations are tested separately beacuse they are the base
     cases of the Bonnet’s recursion algorithm used by mpfr_legendre */
  test_first_iteration ();
  test_second_iteration ();

  /* res precision si arbitrarily low. ARBITRARILY_LOW_PREC should be lower
     than any other precision */
  test_sample_with_precision (IEEE754_SINGLE_PREC, ARBITRARILY_LOW_PREC);
  test_sample_with_precision (IEEE754_DOUBLE_PREC, ARBITRARILY_LOW_PREC);
  test_sample_with_precision (MPFR_PREC_100, ARBITRARILY_LOW_PREC);
  test_sample_with_precision (MPFR_PREC_200, ARBITRARILY_LOW_PREC);

  /* res precision is IEEE754_DOUBLE_PREC. It's higher than
     IEEE754_SINGLE_PREC and lower than MPFR_PREC_100 and
     MPFR_PREC_200 */
  test_sample_with_precision (IEEE754_SINGLE_PREC, IEEE754_DOUBLE_PREC);
  test_sample_with_precision (IEEE754_DOUBLE_PREC, IEEE754_DOUBLE_PREC);
  test_sample_with_precision (MPFR_PREC_100, IEEE754_DOUBLE_PREC);
  test_sample_with_precision (MPFR_PREC_200, IEEE754_DOUBLE_PREC);
  
  /* res precision is MPFR_PREC_200, the highest available.
     NOTE: to reach higher precisions, the expected values should be
     recalculated, because they all contains a 200-bit significand */
  test_sample_with_precision (IEEE754_SINGLE_PREC, MPFR_PREC_200);
  test_sample_with_precision (IEEE754_DOUBLE_PREC, MPFR_PREC_200);
  test_sample_with_precision (MPFR_PREC_100, MPFR_PREC_200);
  test_sample_with_precision (MPFR_PREC_200, MPFR_PREC_200);

  /* This procedure test the rouding of the Legendre algorithm with
     different precisions. There's no need to include special cases
     in the following tests, because they tests their own rounding */
  test_round (IEEE754_SINGLE_PREC);
  test_round (IEEE754_DOUBLE_PREC);
  test_round (MPFR_PREC_100);
  test_round (MPFR_PREC_200);

  tests_end_mpfr ();

  return 0;
}
