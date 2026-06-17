/* Test file for lambertW.

Copyright 2026 Free Software Foundation, Inc.
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

#define INVE "0.367879441171442321595523770161"

static void
w0_test_domain (void)
{
  int ternary;
  mpfr_t inve, res;

  mpfr_init2 (res, IEEE754_DOUBLE_PREC);
  mpfr_init2 (inve, IEEE754_DOUBLE_PREC);
  mpfr_set_str (inve, INVE, 10, MPFR_RNDN);
  MPFR_SET_NEG (inve);

  ternary = mpfr_lambertw0 (res, inve, MPFR_RNDN);
  if (ternary != 0 || mpfr_cmp (res, __gmpfr_mone) != 0)
    {
      printf ("error: W_0(-1/e) = -1, got W_0(-1/e) = ");
      mpfr_out_str (stdout, 10, 0, res, MPFR_RNDN);
      putchar ('\n');
      exit (1);
    }

  mpfr_clear (inve);
  mpfr_clear (res);
}

static void
w1_test_domain (void)
{
  int ternary;
  mpfr_t inve, res;

  mpfr_init2 (res, IEEE754_DOUBLE_PREC);
  mpfr_init2 (inve, IEEE754_DOUBLE_PREC);
  mpfr_set_str (inve, INVE, 10, MPFR_RNDN);
  MPFR_SET_NEG (inve);

  ternary = mpfr_lambertw0 (res, inve, MPFR_RNDN);
  if (ternary != 0 || mpfr_cmp (res, __gmpfr_mone) != 0)
    {
      printf ("error: W_0(-1/e) = -1, got W_0(-1/e) = ");
      mpfr_out_str (stdout, 10, 0, res, MPFR_RNDN);
      putchar ('\n');
      exit (1);
    }

  mpfr_clear (inve);
  mpfr_clear (res);
}

static void
w0_special_cases (void)
{
  int ternary;
  mpfr_t inf, res;

  mpfr_init2 (inf, IEEE754_DOUBLE_PREC);
  mpfr_init2 (res, IEEE754_DOUBLE_PREC);

  mpfr_set_inf (inf, 1);
  ternary = mpfr_lambertw0 (res, inf, MPFR_RNDN);
  if (ternary != 0 || !MPFR_IS_INF (res) || !MPFR_IS_POS (res))
    {
      printf ("error: W_0(+Inf) = +Inf, got W_0(+Inf) = ");
      mpfr_out_str (stdout, 10, 0, res, MPFR_RNDN);
      putchar ('\n');
      exit (1);
    }

  mpfr_clear (inf);
  mpfr_clear (res);
}

static void
w1_special_cases (void)
{
  int ternary;
  mpfr_t zero, res;

  mpfr_init2 (zero, IEEE754_DOUBLE_PREC);
  mpfr_init2 (res, IEEE754_DOUBLE_PREC);

  /* W_{-1}(+0) = -Inf */
  mpfr_set_zero (zero, 1);
  ternary = mpfr_lambertw1 (res, zero, MPFR_RNDN);
  if (ternary != 0 || !MPFR_IS_INF (res) || !MPFR_IS_NEG (res))
    {
      printf ("error: W_{-1}(+0) = -Inf, got W_{-1}(+0) = ");
      mpfr_out_str (stdout, 10, 0, res, MPFR_RNDN);
      putchar ('\n');
      exit (1);
    }

  /* W_{-1}(-0) = -Inf */
  mpfr_set_zero (zero, -1);
  ternary = mpfr_lambertw1 (res, zero, MPFR_RNDN);
  if (ternary != 0 || !MPFR_IS_INF (res) || !MPFR_IS_NEG (res))
    {
      printf ("error: W_{-1}(-0) = -Inf, got W_{-1}(-0) = ");
      mpfr_out_str (stdout, 10, 0, res, MPFR_RNDN);
      putchar ('\n');
      exit (1);
    }

  mpfr_clear (zero);
  mpfr_clear (res);
}

int
main (void)
{
  tests_start_mpfr ();

  /* W_0 */

  /* the domain of W_0 is [-1/e, +Inf) */
  w0_test_domain ();
  w0_special_cases ();

  /* W_{-1} */

  /* the domain of W_0 is [-1/e, 0) */
  w1_test_domain ();
  w1_special_cases ();

  tests_end_mpfr ();
  return 0;
}
