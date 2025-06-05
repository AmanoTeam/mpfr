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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

int
mpfr_legendre (mpfr_ptr res, int n, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int ret = 0;
  unsigned is_within_domain = 1;
  mpfr_t one, minus_one;

  /* the following variables are used (and consequently initialized) only
     for n >= 2, where x is not equal to -1, 0 or 1 */
  unsigned i, basic_ops_bits;
  mpfr_t p1, p2, pn, first_term, second_term;
  mpfr_prec_t res_prec, realprec, test_prec, err;
  MPFR_GROUP_DECL (group);
  MPFR_ZIV_DECL (loop);

  mpfr_init2 (one, MPFR_PREC_MIN);
  mpfr_init2 (minus_one, MPFR_PREC_MIN);
  MPFR_SET_ONE (one);
  MPFR_SET_ONE (minus_one);
  MPFR_SET_NEG(minus_one);

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", MPFR_PREC (x), mpfr_log_prec, x, rnd_mode),
     ("legendre[%Pd]=%.*Rg", MPFR_PREC (res), mpfr_log_prec, res));

  /* first, check if x belongs to the domain [-1,1]. If it's not,
     res is set to NAN, and 0 is returned */
  is_within_domain &= mpfr_lessequal_p(x, one);
  is_within_domain &= mpfr_greaterequal_p (x, minus_one);

  if (!is_within_domain)
    {
      MPFR_SET_NAN (res);
      /* As specified in the documentation, "[...] a NaN result
         (Not-a-Number) always corresponds to an exact return value." */
      ret = 0;
      goto cleanup;
    }

  /* 1 and -1 are the (respectively) upper and lower bound of the Legendre
     polynomial's canonical domain. We can evaluate Pn (for any n) without
     using Bonnet's recursion. Both 1 and -1 are exactly representable,
     so we this will always return 0 */
  if (mpfr_equal_p (x, one))
    {
      mpfr_set_ui (res, 1, rnd_mode);
      ret = 0;
      goto cleanup;
    }
  if (mpfr_equal_p (x, minus_one))
    {
      mpfr_set_si (res, (n&1) == 0 ? 1 : -1, rnd_mode);
      ret = 0;
      goto cleanup;
    }

  /* P_0 = 1 */
  if (n == 0)
    {
      mpfr_set_ui (res, 1, rnd_mode);
      /* 1 is exactly representable in MPFR regardless of precision,
         so this will always return 0 */
      ret = 0;
      goto cleanup;
    }
  /* P_1 = x */
  if (n == 1)
    {
      /* result is set to x. The ternary value of mpfr_set is returned */
      ret = mpfr_set(res, x, rnd_mode);
      goto cleanup;
    }

  /* Pn(0) = 0 if n is odd */
  if (MPFR_IS_ZERO (x) && (n&1) == 1)
    {
      MPFR_SET_ZERO (res);
      /* 0 is exactly representable in MPFR regardless of precision,
         so this will always return 0 */
      ret = 0;
      goto cleanup;
    }

  /* each iteration of the recurrence performs 5 rounded operations:
       - 2 integer multiplications (mpfr_mul_ui),
       - 1 fused multiply-subtract (mpfr_fms),
       - 1 integer division (mpfr_div_ui),
       - 1 assignment with rounding (mpfr_set).
     Hence the per-step multiplicative factor is of the form (1 + theta)^5
     with |theta| <= 2^(-realprec). Expanding the binomial:
     (1 + theta)^5 = 1 + 5theta + 10theta^2 + 10theta^3 + 5theta^4 + theta^5
     Since |theta| <= 2^(-realprec), for realprec >=2, theta <= 1/4.
     5 + 10*1/4 + 10*1/16 + 5*1/64 + 1/256 + 1 = 9.20703125 < 10.
     Since |(1 + theta)^5| <= 10, the error amplification
     of each step is bounded by 10 * |theta|. Therefore, after n steps, the
     total relative error is bounded by:
     total_error <= 2^(n * log2(10) - realprec). For a better stability,
     we rewrite it as total_error <= 2^(n * L - B - 5 - realprec),
     where L = ceil(log2(10)), and B = ceil(log2(5n)). */

  res_prec = MPFR_PREC(res);
  basic_ops_bits = MPFR_INT_CEIL_LOG2(n * 5);
  /* working precision = target + small margin (5 bits) + ceil(log2(10)) ~ 3 */
  realprec = res_prec + basic_ops_bits + 8;
  /* error = ceil(log2(n)) + guard bits + small margin. ceil(log2(10)) ~ 3
     has been added directly to the margin of 5 bits. */
  err = MPFR_INT_CEIL_LOG2(n) + basic_ops_bits + 8;

  MPFR_GROUP_INIT_5 (group, realprec,
                     p1, p2, pn, first_term, second_term);

  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
      i = 2;
      /* p1 = P_1 = x */
      mpfr_set (p1, x, MPFR_RNDN);
      /* p2 = P_0 = 1 */
      mpfr_set_ui (p2, 1, MPFR_RNDN);
      while (i <= n)
        {
          /* below, for each theta<i> holds the following inequality:
             theta<i>| <= 2^(-realprec) */

          mpfr_mul_ui (first_term, x, 2 * i - 1, MPFR_RNDN);
          /* first_term = (x*(2i-1)) * (1 + theta) */
          mpfr_mul_ui (second_term, p2, i - 1, MPFR_RNDN);
          /* second_term = (p2·(i-1)) * (1 + theta2) */
          mpfr_fms (pn, first_term, p1, second_term, MPFR_RNDN);
          /* using fms instead of multiplication + subtraction to perform just
             one rounding for the final result. So,
             pn = (first_term * p1 - second_term) * (1 + theta3) */
          mpfr_div_ui (pn, pn, i, MPFR_RNDN);
          /* pn = (pn_tmp / i) * (1 + theta4) */

          mpfr_swap (p2, p1);
          mpfr_set (p1, pn, MPFR_RNDN);
          /* this set operation introduces a new rounding error bounded by
             2^(-realprec) */
          i++;
        }

      test_prec = realprec - err;

      if (mpfr_min_prec (pn) < test_prec)
        break;

      if (MPFR_CAN_ROUND (pn, test_prec, res_prec, rnd_mode))
        break;

      MPFR_ZIV_NEXT (loop, realprec);
      MPFR_GROUP_REPREC_5 (group, realprec,
                           p1, p2, pn, first_term, second_term);
    }

  MPFR_ZIV_FREE (loop);

  ret = mpfr_set (res, pn, rnd_mode);

  MPFR_GROUP_CLEAR (group);

cleanup:
  mpfr_clear (one);
  mpfr_clear (minus_one);
  return ret;
}
