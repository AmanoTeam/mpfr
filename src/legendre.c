/* legendre -- Compute the nth degree Legendre polynomial.

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

#include "mpfr-impl.h"

int
mpfr_legendre (mpfr_ptr res, int n, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int ret = 0, is_within_domain = 1;
  unsigned i = 2;
  mpfr_prec_t realprec = MPFR_GET_PREC (res) + 20;
  mpfr_t one, minus_one;
  mpfr_init2 (one, MPFR_PREC_MIN);
  mpfr_init2 (minus_one, MPFR_PREC_MIN);
  MPFR_SET_ONE (one);
  MPFR_SET_ONE (minus_one);
  MPFR_SET_NEG(minus_one);

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", mpfr_get_prec (x), mpfr_log_prec, x, rnd_mode),
     ("legendre[%Pd]=%.*Rg ret=%d",
      mpfr_get_prec (res), mpfr_log_prec, res, ret));

  /* First, check if x belongs to the domain [-1,1]. If it's not,
     res is set to NAN */
  is_within_domain &= mpfr_lessequal_p(x, one);
  is_within_domain &= mpfr_greaterequal_p (x, minus_one);

  if (!is_within_domain)
    {
      MPFR_SET_NAN (res);
      goto cleanup;
    }

  /* 1 and -1 are the (respectively) upper and lower bound of the Legendre
     polynomial's canonical domain. We can evaluate Pn (for any n) without
     using Bonnet's recursion */
  if (mpfr_equal_p (x, one))
    {
      mpfr_set_ui (res, 1, rnd_mode);
      goto cleanup;
    }
  if (mpfr_equal_p (x, minus_one))
    {
      mpfr_set_si (res, (n&1) == 0 ? 1 : -1, rnd_mode);
      goto cleanup;
    }

  /* P_0 = 1 */
  if (n == 0)
    {
      mpfr_set_ui (res, 1, rnd_mode);
      goto cleanup;
    }
  /* P_1 = x */
  if (n == 1)
    {
      mpfr_set(res, x, rnd_mode);
      goto cleanup;
    }

  /* Pn(0) = 0 if n is odd */
  if (MPFR_IS_ZERO (x) && (n&1) == 1)
    {
      MPFR_SET_ZERO (res);
      goto cleanup;
    }

  mpfr_t p2;
  mpfr_t p1;
  mpfr_t pn;
  mpfr_t first_term;
  mpfr_t second_term;

  mpfr_init2 (first_term, realprec);
  mpfr_init2 (second_term, realprec);
  mpfr_init2 (p2, realprec);
  mpfr_init2 (p1, realprec);
  mpfr_init2 (pn, realprec);

  mpfr_set (p2, one, MPFR_RNDD);
  mpfr_set (p1, x, MPFR_RNDD);

  while (i <= n)
    {
      /* first_term = x*(2i-1) */
      mpfr_mul_ui (first_term, x, 2*i-1, MPFR_RNDD);
      /* second_term = p2*(i-1) */
      mpfr_mul_ui (second_term, p2, i-1, MPFR_RNDD);
      /* pn = first_term*p1 - second_term */
      mpfr_fms (pn, first_term, p1, second_term, MPFR_RNDD);
      /* pn = pn/i */
      mpfr_div_ui (pn, pn, i, MPFR_RNDD);

      mpfr_swap (p2, p1);
      mpfr_set (p1, pn, MPFR_RNDD);
      i++;
    }

  /* Rounding before mpfr_set guarantee that the final value in res is exact,
     since prec(pn) >> prec(res) */
  mpfr_prec_round (pn, MPFR_GET_PREC (res), rnd_mode);
  ret = mpfr_set (res, pn, rnd_mode);

  mpfr_clear (first_term);
  mpfr_clear (second_term);
  mpfr_clear (p2);
  mpfr_clear (p1);
  mpfr_clear (pn);

cleanup:
  mpfr_clear (one);
  mpfr_clear (minus_one);
  return ret;
}
