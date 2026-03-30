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
  int ternary_value = 0, inex;

  /* the following variables are used (and consequently initialized) only
     for n >= 2, where x is not equal to -1, 0 or 1 */
  long i;
  mpfr_t p1, p2, pn, first_term, second_term;
  mpfr_prec_t res_prec, realprec, x_prec;
  mpfr_exp_t lost_bits;
  mpfr_exp_t b_i, f_i, g_i, h_i, q_i, a_i;
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

  /* H_n(0) when n is an odd number is always 0 */
  if (MPFR_IS_ZERO (x) && (n & 1))
    {
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

  realprec = x_prec > res_prec ? x_prec : res_prec + 10;
  realprec += MPFR_INT_CEIL_LOG2 (realprec);

  MPFR_GROUP_INIT_5 (group, realprec,
                     p1, p2, pn, first_term, second_term);

  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
      i = 1;
      
      /* p1 = 2x, p2 = 1 */
      inex = mpfr_mul_ui (p1, x, 2, MPFR_RNDN);
      mpfr_set_ui (p2, 1, MPFR_RNDN);

      b_i = LONG_MIN;                         /* 2^b_i is the absolute error on p2 */
      a_i = MPFR_GET_EXP (p1) - realprec - 1; /* 2^a_i is the absolute error on p1 */

      while (i < n)
        {
          /* first_term = 2x */
          inex |= mpfr_mul_ui (first_term, x, 2, MPFR_RNDN);
          f_i = MPFR_GET_EXP (first_term) - realprec - 1;

          /* second_term = p2 * 2i */
          inex |= mpfr_mul_ui (second_term, p2, 2 * i, MPFR_RNDN);
          g_i = MAX (MPFR_GET_EXP (second_term) - realprec,
                     b_i + MPFR_INT_CEIL_LOG2 (2*i) + 1);

          /* first_term = first_term * p1 */
          inex |= mpfr_mul (first_term, first_term, p1, MPFR_RNDN);
          h_i = 2 + MAX3 (MPFR_GET_EXP (first_term) - realprec - 1,
                          f_i + MPFR_GET_EXP (p1), 2 + MPFR_GET_EXP (x) + a_i);

          /* pn = first_term - second_term */
          inex |= mpfr_sub (pn, first_term, second_term, MPFR_RNDN);
          q_i = 2 + MAX3 (MPFR_GET_EXP (pn) - realprec - 1, h_i, g_i);

          /* p2 = p1, p1 = pn */
          mpfr_swap (p2, p1);
          mpfr_swap (p1, pn);
          b_i = a_i;
          a_i = q_i;

          i++;
        }

      lost_bits = a_i - (MPFR_GET_EXP (p1) - realprec);

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
