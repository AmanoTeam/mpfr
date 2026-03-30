/* Test file for (physicist's) hermite polynomials.

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

#define MPFR_ORTHOGONAL_POLY_FN mpfr_hermite
#define X_LOWER_BOUND -128.0
#define X_HIGHER_BOUND 128.0
#include "torthopoly.c"

#define RANDOM_TESTS_N_DEGREE 10
#define RANDOM_TESTS_BATCH    20
#define DYADIC_BOUND          10

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

  /* H_3(3.49376) = 2.9924358881463502e2 with MPFR_RNDN */
  mpfr_set_d (x, 3.49376, MPFR_RNDN);
  mpfr_set_str (expected, "299.24358881463500800000000000000", 10, MPFR_RNDN);

  mpfr_hermite (res, 3, x, MPFR_RNDN);
  if (mpfr_cmp (res, expected))
    {
      printf ("test_double_precision failed;\n");
      DUMP_NUMBERS (expected, res);
      exit (1);
    }

  /* H_6(-6.25) = 3.1102803906250000e6 with MPFR_RNDN */
  mpfr_set_d (x, -6.25, MPFR_RNDN);
  mpfr_set_str (expected, "3.11028039062500000000000000000e6", 10, MPFR_RNDN);

  mpfr_hermite (res, 6, x, MPFR_RNDN);
  if (mpfr_cmp (res, expected))
    {
      printf ("test_double_precision failed;\n");
      DUMP_NUMBERS (expected, res);
      exit (1);
    }

  /* H_2(0.0001) = -1.9999999600000000e0 with MPFR_RNDN */
  mpfr_set_d (x, 0.0001, MPFR_RNDN);
  mpfr_set_str (expected, "-1.999999960000000000000000000000", 10, MPFR_RNDN);

  mpfr_hermite (res, 2, x, MPFR_RNDN);
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

/* Exhaustive test of degree-n Hermite polynomial with all fractions
   a/2^b with |a| <= A and 0 <= b <= B, for precision p.
   Assume n >= 1. */
static void
test_exact (int n, int A, int B, mpfr_prec_t p)
{
  mpq_t *H0, *H1, t, u;
  int i, j, a, b, rnd;
  mpfr_t x, y, z;

  H0 = (mpq_t*) malloc ((n + 1) * sizeof (mpq_t));
  H1 = (mpq_t*) malloc ((n + 1) * sizeof (mpq_t));
  for (i = 0; i <= n; i++) {
    mpq_init (H0[i]); /* set to 0 */
    mpq_init (H1[i]); /* set to 0 */
  }
  mpq_init (t);
  mpq_init (u);
  /* use the physicist's Hermite recurrence:
     H_0(x) = 1, H_1(x) = 2x,
     H_j(x) = 2x * H_{j-1}(x) - 2(j-1) * H_{j-2}(x)
     In coefficient form:
     H_j[i] = 2 * H_{j-1}[i-1] - 2*(j-1) * H_{j-2}[i] */
  mpq_set_ui (H0[0], 1, 1); /* H_0 = 1 */
  mpq_set_ui (H1[1], 2, 1); /* H_1 = 2x */
  for (j = 2; j <= n; j++) {
    /* H[j] = 2x * H[j-1] - 2(j-1) * H[j-2]
       thus H[j][i] = 2 * H[j-1][i-1] - 2*(j-1) * H[j-2][i].
       Invariant: H[j-2] is stored in H0, and H[j-1] in H1. */
    for (i = 0; i <= j; i++) {
      if (i == 0)
        mpq_set_ui (t, 0, 1);
      else {
        mpq_set_ui (t, 2, 1);
        mpq_mul (t, t, H1[i-1]);
      }
      /* t = 2 * H[j-1][i-1] */
      mpq_set_ui (u, 2*(j-1), 1);
      mpq_mul (u, u, H0[i]);
      /* u = 2*(j-1) * H[j-2][i] */
      mpq_sub (H0[i], t, u);
      /* now H0[i] contains H[j][i] */
    }
    /* swap H0 and H1 */
    for (i = 0; i <= j; i++)
      mpq_swap (H0[i], H1[i]);
  }

  mpfr_init2 (x, 64);
  mpfr_init2 (y, p);
  mpfr_init2 (z, p);

  for (a = -A; a <= A; a++)
    for (b = 0; b <= B; b++) {
      /* compute t = Hn(a/2^b) */
      mpq_set_si (u, a, 1ul<<b);
      mpq_set (t, H1[n]);
      for (i = n-1; i >= 0; i--) {
        mpq_mul (t, t, u);
        mpq_add (t, t, H1[i]);
      }

      /* now t = Hn(a/2^b) exactly */

      mpfr_set_si_2exp (x, a, -b, MPFR_RNDN);
      RND_LOOP (rnd) {
        mpfr_rnd_t r = (mpfr_rnd_t) rnd;
        mpfr_set_q (y, t, (mpfr_rnd_t) rnd); /* expected result */
        mpfr_hermite (z, n, x, r);
        if (mpfr_cmp (y, z)) {
          printf ("Error in test_exact for n=%d a=%d b=%d p=%lu rnd=%s\n",
                  n, a, b, p, mpfr_print_rnd_mode (r));
          DUMP_NUMBERS (y, z);
          exit (1);
        }
      }
    }

  for (i = 0; i <= n; i++) {
    mpq_clear (H0[i]);
    mpq_clear (H1[i]);
  }
  free (H0);
  free (H1);
  mpq_clear (t);
  mpq_clear (u);
  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

static void
test_exact_dyadic (void)
{
  int n;
  mpfr_prec_t p;

  for (n = 1; n <= 10; n++)
    for (p = DYADIC_BOUND - 3; p <= DYADIC_BOUND; p++)
      test_exact (n, DYADIC_BOUND, DYADIC_BOUND, p);
}

int
main (void)
{
  tests_start_mpfr ();

  /* for singular input (+/-Inf or NaN), mpfr_hermit should always set res
     to NaN and return 0 */
  test_singular_input ();

  /* the first two iterations are tested separately because they are the base
     cases of the recursion algorithm used to calculate mpfr_hermite */
  test_first_iteration ();
  test_second_iteration ();

  /* test physicist's Hermite polynomials with fixed test cases in double
     precision */
  test_double_precision ();

  random_poly_suite (RANDOM_TESTS_N_DEGREE, RANDOM_TESTS_BATCH,
                     IEEE754_DOUBLE_PREC);

  test_exact_dyadic ();

  tests_end_mpfr ();
  return 0;
}

#undef MPFR_ORTHOGONAL_POLY_FN
#undef X_LOWER_BOUND
#undef X_HIGHER_BOUND
#undef RANDOM_TESTS_N_DEGREE
#undef RANDOM_TESTS_BATCH
