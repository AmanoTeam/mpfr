/* hermite -- Compute the nth degree (physicist's) Hermite polynomial.

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

int
mpfr_hermite (mpfr_ptr res, unsigned n, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int ternary_value = 0;

  /* the following variables are used (and consequently initialized) only
     for n >= 2, where x is not equal to -1, 0 or 1 */
  // unsigned i;
  // mpfr_t p1, p2, pn, first_term, second_term;
  // mpfr_prec_t res_prec, realprec, x_prec, test_prec, err;
  // MPFR_GROUP_DECL (group);
  // MPFR_ZIV_DECL (loop);
  mpfr_prec_t res_prec, realprec;

  res_prec = MPFR_PREC (res);

  /* */
  if (n > 8192)
    goto nan_ret;

  /* H_0(x) = 1. In this case, since the output is const and does not depend
     on the value of x, no further analysis on the value of x is performed */
  if (n == 0)
    {
      mpfr_set_ui (res, 1, rnd_mode);
      /* 1 is exactly representable in MPFR regardless of precision,
        so this will always return 0 */
      return 0;
    }

  if (MPFR_IS_NAN (x) || MPFR_IS_INF (x))
    {
  nan_ret:
      MPFR_SET_NAN (res);
      /* as specified in the documentation, "[...] a NaN result
         (Not-a-Number) always corresponds to an exact return value." */
      return 0;
    }

  if (MPFR_IS_ZERO (x))
    {
      /* H_n(0) when n is an even number (i.e. n=2k) is -1^k * (2k!)/k!.
         Here we are computing it with the Log-Gamma method to avoid overflow
         for large n, so we have exp(lngamma(2*k + 1) - lngamma(k + 1)), with
         k = n / 2. see algorithms.tex for further details */
      if ((n&1) == 0)
        {
          mpfr_t first, second, gamma1, gamma2, sub, e;
          MPFR_GROUP_DECL (group);
          MPFR_BLOCK_DECL (flags);

          realprec = res_prec + 10;

          MPFR_GROUP_INIT_6 (group, realprec, first, second, gamma1, gamma2,
                             sub, e);

          mpfr_set_ui (first, n + 1, MPFR_RNDN);
          mpfr_set_ui (second, (n >> 1) + 1, MPFR_RNDN);
          mpfr_lngamma (gamma1, first, MPFR_RNDN);
          mpfr_lngamma (gamma2, second, MPFR_RNDN);
          mpfr_sub (sub, gamma1, gamma2, MPFR_RNDN);

          /* exp can overflow */
          MPFR_BLOCK (flags, mpfr_exp (e, sub, MPFR_RNDN));

          /* in case of overflow, res is set to NaN, and 0 is returned */
          if (MPFR_OVERFLOW (flags))
            goto nan_ret;

          ternary_value = mpfr_set (res, e, rnd_mode);

          /* -1^k is negative iff k = n/2 is odd */
          if (((n >> 1) & 1))
            MPFR_SET_NEG (res);

          MPFR_GROUP_CLEAR (group);

          return ternary_value;
        }

        /* H_n(0) when n is an odd number is always 0 */
        if (n&1)
          {
            MPFR_SET_ZERO (res);
            /* 0 is exactly representable in MPFR regardless of precision,
               so this will always return 0 */
            return 0;
          }
    }

  /* H_1(x) = 2x */
  if (n == 1)
    {
      /* result is set to 2x. The ternary value of mpfr_set is returned */
      return mpfr_mul_ui (res, x, 2, rnd_mode);
    }

  return ternary_value;
}
