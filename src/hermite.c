/* hermite -- Compute the nth degree (physicist's) Hermite polynomial.

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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

static int
zero_x_even_degree (mpfr_ptr res, unsigned n, mpfr_rnd_t rnd_mode)
{
  int ternary_value;
  mpfr_exp_t exp_sub;
  mpfr_t first, second, gamma1, gamma2, sub, e;
  mpfr_prec_t res_prec, realprec, test_prec, err;
  MPFR_ZIV_DECL (loop);
  MPFR_GROUP_DECL (group);
  MPFR_BLOCK_DECL (flags);

  res_prec = MPFR_PREC (res);
  realprec = res_prec + 10;
  /* we start */
  err = (mpfr_exp_t) MPFR_INT_CEIL_LOG2 (n) + 2;
  MPFR_GROUP_INIT_6 (group, realprec, first, second, gamma1,
                      gamma2, sub, e);

  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
      mpfr_set_ui (first, n + 1, MPFR_RNDN);
      mpfr_set_ui (second, (n >> 1) + 1, MPFR_RNDN);
      mpfr_lngamma (gamma1, first, MPFR_RNDN);
      mpfr_lngamma (gamma2, second, MPFR_RNDN);
      mpfr_sub (sub, gamma1, gamma2, MPFR_RNDN);

      /* exp may overflow */
      MPFR_BLOCK (flags, mpfr_exp (e, sub, MPFR_RNDN));

      /* in case of overflow, res is set to NaN, and 0 is returned. We skip
         both the last set and the sign set by jumping to mpfr_clear */
      if (MPFR_OVERFLOW (flags))
        {
          MPFR_SET_NAN (res);
          ternary_value = 0;
          goto clear;
        }

      /* in algorithms.tex there is a detail error analysis for this algorithm.
         Here we use a practical heuristic like the following.
         Let s be the exact value and sub = s + ds its MPFR approximation,
         with |ds| <= k3 * ulp(s) = k3 * 2^(Exp(s) - w), where k3 = 3/2.
         For w large enough, we have |ds| <= 1/2, hence
            |exp(ds) - 1| <= 4/3 * |ds|.
         Therefore,
             exp(sub) = exp(s) * (1 + 2^(Exp(s)+2-w)) up to rounding.
         Including the final rounding error of mpfr_exp, we get
             e = exp(s) * (1 + 2^(err-w)),
         with:
             err = Exp(s) + 3   if Exp(s)+2 > 0;
                   2            if Exp(s)+2 = 0;
                   1            if Exp(s)+2 < 0. */
      exp_sub = MPFR_GET_EXP (sub);
      err = (exp_sub + 2 > 0)
            ? exp_sub + 3
            : (exp_sub + 2 == 0) ? 2 : 1;
      test_prec = realprec - err;

      if (mpfr_min_prec (e) < test_prec - 1)
        break;

      if (MPFR_LIKELY (MPFR_CAN_ROUND (e, test_prec, res_prec, rnd_mode)))
        break;
    }

  ternary_value = mpfr_set (res, e, rnd_mode);

  /* -1^k is negative iff k = n/2 is odd */
  if (((n >> 1) & 1))
    MPFR_SET_NEG (res);

clear:
  MPFR_ZIV_FREE (loop);
  MPFR_GROUP_CLEAR (group);

  return ternary_value;
}

static mpfr_prec_t
error_extra_bits (unsigned n)
{
  /* pre-computed error bits for n <= 10, see the table in algorithms.tex */
  const mpfr_prec_t extra_bits[] = { 0, 0, 3, 4, 6, 7, 9, 10, 12, 13, 14 };
  double log2_r1, log2_A, asym_bits;

  if (n <= 10)
    return extra_bits[n];

  /* for larger n we use the asymptotic analysis shown in algorithms.tex,
     in particular we have:
       - log2_r1 = log2(1*sqrt(3)) = 1.4499843134764958;
       - log2_A = log2((4+sqrt(3))/sqrt(3)) = 1.726570147010381;
       - asym_bits = log2_A + n*(log2_r1) + 1 */
  log2_r1 = 1.4499843134764958;
  log2_A = 1.726570147010381;
  asym_bits = log2_A + n * log2_r1 + 1.0;

  return MPFR_CEIL (asym_bits);
}

int
mpfr_hermite (mpfr_ptr res, long n, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int ternary_value = 0;

  /* the following variables are used (and consequently initialized) only
     for n >= 2, where x is not equal to -1, 0 or 1 */
  unsigned i;
  mpfr_t p1, p2, pn, first_term, second_term;
  mpfr_prec_t res_prec, realprec, x_prec, test_prec, err;
  MPFR_GROUP_DECL (group);
  MPFR_ZIV_DECL (loop);

  x_prec = MPFR_PREC (x);
  res_prec = MPFR_PREC (res);

  /* NaN are checke *before* any other check, according to C++ specs:
     "If the argument is NaN, NaN is returned [...]".
     We extended to +/-Inf as well. */
  if (MPFR_IS_NAN (x) || MPFR_IS_INF (x))
    {
      MPFR_SET_NAN (res);
      /* as specified in the documentation, "[...] a NaN result
         (Not-a-Number) always corresponds to an exact return value." */
      return 0;
    }

  /* H_0(x) = 1. In this case, since the output is const and does not depend
     on the value of x, no further analysis on the value of x is performed */
  if (n == 0)
    {
      mpfr_set_ui (res, 1, rnd_mode);
      /* 1 is exactly representable in MPFR regardless of precision,
        so this will always return 0 */
      return 0;
    }

  if (MPFR_IS_ZERO (x))
    {
      /* H_n(0) when n is an even number (i.e. n=2k) is -1^k * (2k!)/k!.
         Here we are computing it with the Log-Gamma method to avoid overflow
         for large n, so we have exp(lngamma(2*k + 1) - lngamma(k + 1)), with
         k = n / 2. see algorithms.tex for further details */
      if ((n&1) == 0)
        return zero_x_even_degree (res, n, rnd_mode);

      /* H_n(0) when n is an odd number is always 0 */
      MPFR_SET_ZERO (res);
      /* 0 is exactly representable in MPFR regardless of precision,
          so this will always return 0 */
      return 0;
    }

  /* H_1(x) = 2x */
  if (n == 1)
    {
      /* result is set to 2x. The ternary value of mpfr_set is returned */
      return mpfr_mul_ui (res, x, 2, rnd_mode);
    }

  /* if x_prec > res_prec, then we use x_prec as the starting precision
     for the Ziv's loop, otherwise we use res_prec. We add then a coefficient
     that is proportional either to x_prec or res_prec + 10 safety bits */
  realprec = x_prec > res_prec ? x_prec : res_prec + 10;
  realprec += MPFR_INT_CEIL_LOG2 (realprec);

  /* the error is set to error_extra_bits, plus 2 extra bits of safety for
     the final rounding */
  err = error_extra_bits (n) + 2;

  MPFR_GROUP_INIT_5 (group, realprec,
                     p1, p2, pn, first_term, second_term);

  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
      i = 1;

      /* p1 = 2x, p2 = 1 */
      mpfr_mul_ui (p1, x, 2, MPFR_RNDN);
      mpfr_set_ui (p2, 1, MPFR_RNDN);
      while (i < n)
        {
          /* first_term = 2x */
          mpfr_mul_ui (first_term, x, 2, MPFR_RNDN);
          /* first_term *= p1 */
          mpfr_mul (first_term, first_term, p1, MPFR_RNDN);
          /* second_term = p2 * 2i */
          mpfr_mul_ui (second_term, p2, 2 * i, MPFR_RNDN);
          /* pn = first_term - second_term */
          mpfr_sub (pn, first_term, second_term, MPFR_RNDN);

          /* p2 = p1, p1 = pn */
          mpfr_set (p2, p1, MPFR_RNDN);
          mpfr_set (p1, pn, MPFR_RNDN);

          i++;
        }

      test_prec = realprec - err;

      if (mpfr_min_prec (pn) < test_prec - 1)
        break;

      if (MPFR_LIKELY (MPFR_CAN_ROUND (pn, test_prec, res_prec, rnd_mode)))
        break;

      MPFR_ZIV_NEXT (loop, realprec);
      MPFR_GROUP_REPREC_5 (group, realprec,
                           p1, p2, pn, first_term, second_term);
    }
  MPFR_ZIV_FREE (loop);
  ternary_value = mpfr_set (res, pn, rnd_mode);

  MPFR_GROUP_CLEAR (group);

  return ternary_value;
}
