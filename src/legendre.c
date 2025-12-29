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

#include <math.h> /* for ceil function */

static mpfr_prec_t
error_extra_bits (unsigned n)
{
  /* pre-computed error bits for n <= 10, see the table in algorithms.tex */
  const mpfr_prec_t extra_bits[] = { 0, 0, 6, 8, 11, 13, 15, 17, 20, 22, 24 };
  double log2_r1, log2_A, asym_bits;

  if (n <= 10)
    return extra_bits[n];

  /* for larger n we use the asymptotic analysis shown in algorithms.tex,
     in particular we have:
       - log2_r1 = log2(2+2*sqrt(2)) = 2.271553303163612;
       - log2_A = log2((7+5*sqrt(2))/4) = 1.815573220150449;
       - asym_bits = log2_A + n*(log2_r1) + 1 */
  log2_r1 = 2.271553303163612;
  log2_A = 1.815573220150449;
  asym_bits = log2_A + n * log2_r1 + 1.0;

  return ceil (asym_bits);
}

int
mpfr_legendre (mpfr_ptr res, unsigned n, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int ternary_value = 0;
  unsigned is_within_domain = 1;

  /* the following variables are used (and consequently initialized) only
     for n >= 2, where x is not equal to -1, 0 or 1 */
  unsigned i;
  mpfr_t p1, p2, pn, first_term, second_term;
  mpfr_prec_t res_prec, realprec, x_prec, test_prec, err;
  MPFR_GROUP_DECL (group);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", MPFR_PREC (x), mpfr_log_prec, x, rnd_mode),
     ("legendre[%Pd]=%.*Rg", MPFR_PREC (res), mpfr_log_prec, res));

  /* first, check if x belongs to the domain [-1,1] and 0 <= n <= 2^13 = 8192.
     If it's not, res is set to NAN, and 0 is returned */
  is_within_domain &= mpfr_lessequal_p(x, __gmpfr_one);
  is_within_domain &= mpfr_greaterequal_p (x, __gmpfr_mone);

  if (!is_within_domain || n > 8192)
    {
      MPFR_SET_NAN (res);
      /* as specified in the documentation, "[...] a NaN result
         (Not-a-Number) always corresponds to an exact return value." */
      return 0;
    }

  /* 1 and -1 are the (respectively) upper and lower bound of the Legendre
     polynomial's canonical domain. We can evaluate Pn (for any n) without
     using Bonnet's recursion. Both 1 and -1 are exactly representable,
     so we this will always return 0 */
  if (mpfr_equal_p (x, __gmpfr_one))
    {
      mpfr_set_ui (res, 1, rnd_mode);
      return 0;
    }
  if (mpfr_equal_p (x, __gmpfr_mone))
    {
      mpfr_set_si (res, (n&1) == 0 ? 1 : -1, rnd_mode);
      return 0;
    }

  /* P_0 = 1 */
  if (n == 0)
    {
      mpfr_set_ui (res, 1, rnd_mode);
      /* 1 is exactly representable in MPFR regardless of precision,
         so this will always return 0 */
      return 0;
    }
  /* P_1 = x */
  if (n == 1)
    {
      /* result is set to x. The ternary value of mpfr_set is returned */
      return mpfr_set(res, x, rnd_mode);
    }

  /* Pn(0) = 0 if n is odd */
  if (MPFR_IS_ZERO (x) && (n&1) == 1)
    {
      MPFR_SET_ZERO (res);
      /* 0 is exactly representable in MPFR regardless of precision,
         so this will always return 0 */
      return 0;
    }

  res_prec = MPFR_PREC (res);
  x_prec = MPFR_PREC (x);
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
      i = 2;

      /* p1 = x, p2 = 1 */
      mpfr_set (p1, x, MPFR_RNDN);
      mpfr_set_ui (p2, 1, MPFR_RNDN);

      while (i <= n)
        {
          /* first_term = x * (2 * i - 1) */
          mpfr_mul_ui (first_term, x, 2 * i - 1, MPFR_RNDN);
          /* first_term *= p1 */
          mpfr_mul (first_term, first_term, p1, MPFR_RNDN);
          /* second_term = p2 * (i - 1) */
          mpfr_mul_ui (second_term, p2, i - 1, MPFR_RNDN);
          /* pn = first_term - second_term */
          mpfr_sub (pn, first_term, second_term, MPFR_RNDN);
          /* pn = pn/i */
          mpfr_div_ui (pn, pn, i, MPFR_RNDN);

          /* p2 = p1, p1 = pn */
          mpfr_set (p2, p1, MPFR_RNDN);
          mpfr_set (p1, pn, MPFR_RNDN);

          i++;
        }

      /* we have (eventually) added the lost bits due to cancellation
         above, so test_prec is the effective precision we can use */
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
