/* legendre -- Compute the nth degree Legendre polynomial.

Copyright 2025-2026 Free Software Foundation, Inc.
Contributed by Matteo Nicoli and Paul Zimmermann.

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
mpfr_legendre (mpfr_ptr res, long n, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int ternary_value = 0;
  unsigned is_within_domain = 1;

  /* the following variables are used (and consequently initialized) only
     for n >= 2, where x is not equal to -1, 0 or 1 */
  long i;
  mpfr_t p1, p2, pn, first_term, second_term;
  mpfr_prec_t res_prec, realprec;
  mpfr_exp_t lost_bits;
  mpfr_exp_t b_i, log2_i_m1, f_i, g_i, h_i, q_i, a_i, a_n;
  int inex;

  MPFR_GROUP_DECL (group);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", MPFR_PREC (x), mpfr_log_prec, x, rnd_mode),
     ("legendre[%Pd]=%.*Rg", MPFR_PREC (res), mpfr_log_prec, res));

  MPFR_ASSERTN(n >= 0); /* check n is non-negative */

  /* first, check if x belongs to the domain [-1,1].
     If it's not, res is set to NAN, and 0 is returned */
  is_within_domain &= mpfr_lessequal_p (x, __gmpfr_one);
  is_within_domain &= mpfr_greaterequal_p (x, __gmpfr_mone);

  if (!is_within_domain)
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
  /* Analyzing all the test cases where the result is not exact (inex != 0),
     we find that the average number of bits lost per iteration, i.e.,
     lost_bits/(n-1), is about 3.82. We thus add 4*n guard bits.
     For revision 94a5659, we have a total of 615575 such tests.
     With 4n+12 below, we get a probability of failure of 4.7%.
     With 4n+20, we get a probability of failure of 0.7%. */
  realprec = res_prec + 4 * n + 20;
  realprec += MPFR_INT_CEIL_LOG2 (realprec);

  MPFR_GROUP_INIT_5 (group, realprec,
                     p1, p2, pn, first_term, second_term);

  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
      i = 2;

      /* p1 = x, p2 = 1 */
      inex = mpfr_set (p1, x, MPFR_RNDN);     /* p1 is a is algorithms.tex */
      mpfr_set_ui (p2, 1, MPFR_RNDN);         /* exact, p2 is b */
      b_i = LONG_MIN;                         /* 2^b_i is the absolute error on p2 */
      a_i = MPFR_GET_EXP (p1) - realprec - 1; /* 2^a_i is the absolute error on p1 */

      while (i <= n)
        {
          log2_i_m1 = MPFR_INT_CEIL_LOG2(i-1);

          /* first_term = x * (2 * i - 1), with absolute error at step i
             (denoted f_i in algorithms.tex)
             bounded by f_i <= exp(first_term) - p - 1 */
          inex |= mpfr_mul_ui (first_term, x, 2 * i - 1, MPFR_RNDN);
          f_i = MPFR_GET_EXP (first_term) - realprec - 1;

          /* second_term = p2 * (i - 1), with absolute error at step i bounded by
             g_i <= max(exp(second_term)-p, error(p2) + MPFR_INT_CEIL_LOG2(i-1)+1)
           */
          inex |= mpfr_mul_ui (second_term, p2, i - 1, MPFR_RNDN);
          g_i = MAX (MPFR_GET_EXP (second_term) - realprec, b_i + log2_i_m1 + 1);

          /* first_term = first_term * p1, with absolute error at step i
             bounded by
             h_i <= 2 + max(exp(first_term)-p-1, f_i+exp(p1),
                            MPFR_INT_CEIL_LOG2(2*i-1)+MPFR_GET_EXP(x)+a_i) */
          inex |= mpfr_mul (first_term, first_term, p1, MPFR_RNDN);
          h_i = 2 + MAX3 (MPFR_GET_EXP (first_term) - realprec - 1,
                          f_i + MPFR_GET_EXP (p1),
                          MPFR_INT_CEIL_LOG2 (2*i-1) + MPFR_GET_EXP (x) + a_i);

          /* pn = first_term - second_term, with absolute error at step i
             bounded by
             q_i <= 2 + max(exp(pn)-p-1, error(first_term), error(second_term)) */
          inex |= mpfr_sub (pn, first_term, second_term, MPFR_RNDN);
          q_i = 2 + MAX3 (MPFR_GET_EXP (pn) - realprec - 1, h_i, g_i);

          /* pn = pn/i, with absolute error at step i
             bounded by a_i <= max(exp(pn)-p, q_i-MPFR_INT_CEIL_LOG2(i-1)+2) */
          inex |= mpfr_div_ui (pn, pn, i, MPFR_RNDN);
          a_n = MAX (MPFR_GET_EXP (pn) - realprec, q_i - log2_i_m1 + 2);

          /* p2 = p1, p1 = pn */
          mpfr_swap (p2, p1); /* now p2 approximates P_{i-1}(x) */
          mpfr_swap (p1, pn); /* now p1 approximates P_i(x) */
          b_i = a_i;          /* 2^b_i is a bound on the absolute error on p2 */
          a_i = a_n;          /* 2^a_i is a bound on the absolute error on p1 */

          i++;
        }

      /* Now p1 approximates P_n(x), and 2^a_i is a bound on its absolute error.
         Since ulp(p1) = 2^(EXP(p1)-realprec)
         we get the relative error is bounded by:
         2^(a_i - (EXP(p1) - realprec - 1)) */
      lost_bits = a_i - (MPFR_GET_EXP (p1) - realprec);

      /* if inex=0, then all the computation was exact, thus p1 is exactly P_n(x),
         otherwise we call MPFR_CAN_ROUND() to check if we can deduce the correct
         rounding */
      if (inex == 0 || (lost_bits < realprec &&
                     MPFR_CAN_ROUND (p1, realprec - lost_bits, res_prec, rnd_mode)))
        break;
      MPFR_ZIV_NEXT (loop, realprec);
      MPFR_GROUP_REPREC_5 (group, realprec,
                           p1, p2, pn, first_term, second_term);
    }
  MPFR_ZIV_FREE (loop);
  ternary_value = mpfr_set (res, p1, rnd_mode);

  MPFR_GROUP_CLEAR (group);

  return ternary_value;
}
