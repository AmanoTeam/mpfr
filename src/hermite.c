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

/* max (x, y, z) */
#define MAX3(x,y,z) (MAX (x, MAX (y, z)))

int
mpfr_hermite (mpfr_ptr res, long n, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  long i;
  int ternary_value = 0, inex;
  mpfr_t p1, p2, pn, first_term, second_term;
  mpfr_prec_t res_prec, realprec, guard_bits;
  mpfr_exp_t lost_bits;
  mpfr_exp_t b_i, f_i, g_i, h_i, q_i, a_i;
  MPFR_GROUP_DECL (group);
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  res_prec = MPFR_PREC (res);

  /* NaN are checke *before* any other check, according to C++ specs:
     "If the argument is NaN, NaN is returned [...]".
     We extended to +/-Inf as well. */
  if (MPFR_IS_NAN (x) || MPFR_IS_INF (x))
    {
      MPFR_SET_NAN (res);
      /* as specified in the documentation, "[...] a NaN result
         (Not-a-Number) always corresponds to an exact return value." */
      MPFR_RET_NAN;
    }

  /* H_0(x) = 1. In this case, since the output is const and does not depend
     on the value of x, no further analysis on the value of x is performed */
  if (n == 0)
    {
      mpfr_set_ui (res, 1, rnd_mode);
      /* 1 is exactly representable in MPFR regardless of precision,
        so this will always return 0 */
      MPFR_RET (0);
    }

  /* H_n(0) when n is an odd number is always 0 */
  if (MPFR_IS_ZERO (x) && (n & 1))
    {
      MPFR_SET_ZERO (res);
      /* 0 is exactly representable in MPFR regardless of precision,
          so this will always return 0 */
      MPFR_RET (0);
    }

  MPFR_SAVE_EXPO_MARK (expo);

  /* H_1(x) = 2x */
  if (n == 1)
    {
      /* result is set to 2x. The ternary value of mpfr_set is returned */
      ternary_value = mpfr_mul_ui (res, x, 2, rnd_mode);
      goto end;
    }

  /* Analyzing all the test cases where the result is not exact (inex != 0),
     we find that the average number of bits lost per iteration, i.e.,
     lost_bits/(n-1), is about 3.27, but up to about 5.5 for n >= 20.
     We thus add 4*n guard bits for n < 20, and 6*n for n >= 20.
     For revision see eb17cda, where we have a total of 8722 such tests.
     With guard_bits * n + 10, we get a probability of failure of 0.2% */
  guard_bits = n < 20 ? 4 : 6;
  realprec = res_prec + guard_bits * n + 10;
  realprec += MPFR_INT_CEIL_LOG2 (realprec);

  MPFR_GROUP_INIT_5 (group, realprec,
                     p1, p2, pn, first_term, second_term);

  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
      MPFR_BLOCK_DECL (flags);

      i = 1;

      MPFR_BLOCK (flags, inex = mpfr_mul_ui (p1, x, 2, MPFR_RNDN));
      if (MPFR_OVERFLOW (flags))
        {
          /* 2x overflows in extended exponent range;
             H_n(x) overflows for all n >= 1. The sign of H_n(x) for
             large |x| is that of its leading term (2x)^n. */
          ternary_value = mpfr_overflow (res, rnd_mode,
                                         (n & 1) ? MPFR_SIGN (x)
                                                 : MPFR_SIGN_POS);
          MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, MPFR_FLAGS_OVERFLOW);
          break;
        }
      mpfr_set_ui (p2, 1, MPFR_RNDN);           /* exact */

      /* In the loop:
           2^a_i is a bound on the absolute error on p1.
           2^b_i is a bound on the absolute error on p2.
         a_i and b_i come from the previous iteration, and initialized
         below for the first iteration (i = 1). */

      /* 2^b_i is the absolute error on p2 */
      b_i = LONG_MIN;
      /* 2^a_i is the absolute error on p1 */
      a_i = MPFR_GET_EXP (p1) - realprec - 1;

      while (i < n)
        {
          /* first_term = 2x, with absolute error at step i
             (denoted f_i in algorithms.tex)
             bounded by f_i <= exp(first_term) - p - 1 */
          MPFR_BLOCK (flags,
                      inex |= mpfr_mul_ui (first_term, x, 2, MPFR_RNDN));
          if (MPFR_OVERFLOW (flags))
            break;
          f_i = MPFR_GET_EXP (first_term) - realprec - 1;

          /* second_term = p2 * 2i, with absolute error at step i
             bounded by
             g_i <= max(exp(second_term)-p,
                        error(p2) + MPFR_INT_CEIL_LOG2(2*i)+1) */
          MPFR_BLOCK (flags,
                      inex |= mpfr_mul_ui (second_term, p2, 2 * i, MPFR_RNDN));
          if (MPFR_OVERFLOW (flags))
            break;
          g_i = MAX (MPFR_GET_EXP (second_term) - realprec,
                     b_i + MPFR_INT_CEIL_LOG2 (2*i) + 1);

          /* first_term = first_term * p1, with absolute error at step i
             bounded by
             h_i <= 2 + max(exp(first_term)-p-1, f_i+exp(p1),
                            2+MPFR_GET_EXP(x)+a_i) */
          MPFR_BLOCK (flags,
                      inex |= mpfr_mul (first_term, first_term, p1, MPFR_RNDN));
          if (MPFR_OVERFLOW (flags))
            break;
          h_i = 2 + MAX3 (MPFR_GET_EXP (first_term) - realprec - 1,
                          f_i + MPFR_GET_EXP (p1), 2 + MPFR_GET_EXP (x) + a_i);

          /* pn = first_term - second_term, with absolute error at step i
             bounded by
             q_i <= 2 + max(exp(pn)-p-1,
                            error(first_term), error(second_term))
             Note: mpfr_sub can overflow when first_term and second_term
             have opposite signs */
          MPFR_BLOCK (flags,
                      inex |= mpfr_sub (pn, first_term, second_term, MPFR_RNDN));
          if (MPFR_OVERFLOW (flags))
            break;
          q_i = 2 + MAX3 (MPFR_GET_EXP (pn) - realprec - 1, h_i, g_i);

          /* p2 = p1, p1 = pn */
          mpfr_swap (p2, p1); /* now p2 approximates H_{i}(x) */
          mpfr_swap (p1, pn); /* now p1 approximates H_{i+1}(x) */
          b_i = a_i;          /* 2^b_i is a bound on the absolute error on p2 */
          a_i = q_i;          /* 2^a_i is a bound on the absolute error on p1 */

          i++;
        }

      /* If an overflow occurred in the recurrence (detected via flags),
         since we are in extended exponent range, H_n(x) truly overflows.
         The sign is that of the leading term (2x)^n. */
      if (MPFR_OVERFLOW (flags))
        {
          ternary_value = mpfr_overflow (res, rnd_mode,
                                         (n & 1) ? MPFR_SIGN (x)
                                                 : MPFR_SIGN_POS);
          MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, MPFR_FLAGS_OVERFLOW);
          break;
        }

      /* Now p1 approximates H_n(x), and 2^a_i is a bound on its absolute
         error. Since ulp(p1) = 2^(EXP(p1)-realprec), we get the relative
         error is bounded by 2^(a_i - (EXP(p1) - realprec - 1)). */
      lost_bits = a_i - (MPFR_GET_EXP (p1) - realprec);

      /* if inex=0, then all the computation was exact, thus p1 is exactly
         H_n(x), otherwise we call MPFR_CAN_ROUND() to check if we can
         deduce the correct rounding */
      if (inex == 0 ||
          (lost_bits < realprec &&
            MPFR_CAN_ROUND (p1, realprec - lost_bits, res_prec, rnd_mode)))
        {
          ternary_value = mpfr_set (res, p1, rnd_mode);
          break;
        }

      MPFR_ZIV_NEXT (loop, realprec);
      MPFR_GROUP_REPREC_5 (group, realprec,
                           p1, p2, pn, first_term, second_term);
    }
  MPFR_ZIV_FREE (loop);

  MPFR_GROUP_CLEAR (group);

 end:
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (res, ternary_value, rnd_mode);
}
