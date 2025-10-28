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

/* return the number of bits we (eventually) have to add to realprec
   to compute x - y exactly, i.e to avoid catastrophic cancellation.
   Return 0 if cancellation is negligible. */
static mpfr_prec_t
likely_cancellation (mpfr_srcptr x, mpfr_srcptr y, mpfr_prec_t prec)
{
  long lost_bits;
  mpfr_exp_t ex, ey, ed, maxexp;
  mpfr_t diff;

  if (MPFR_IS_ZERO (x) || MPFR_IS_ZERO (y))
    return 0;

  if (MPFR_SIGN (x) != MPFR_SIGN (y))
    return 0;

  ex = MPFR_GET_EXP (x);
  ey = MPFR_GET_EXP (y);

  /* magnitudes too different, subtraction is safe */
  if (labs (ex - ey) > 2)
    return 0;

  /* we compute the difference with a slightly higher precision.
     We add extra 8 bits to the working precision to give some guard space
     for the subtraction itself. */
  mpfr_init2 (diff, prec + 8);
  mpfr_sub (diff, x, y, MPFR_RNDN);

  if (MPFR_IS_ZERO (diff))
    {
      mpfr_clear (diff);
      return prec;
    }

  ed = MPFR_GET_EXP (diff);
  maxexp = (ex > ey) ? ex : ey;
  lost_bits = maxexp - ed;

  mpfr_clear (diff);

  /* if lost_bits <= 0, subtraction doesn’t lose any precision.
     If the subtraction loses only one bit of precision, treat it
     as negligible */
  if (lost_bits <= 1)
    return 0;

  return (mpfr_prec_t) lost_bits;
}

int
mpfr_legendre (mpfr_ptr res, int n, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int ternary_value = 0;
  unsigned is_within_domain = 1;

  /* the following variables are used (and consequently initialized) only
     for n >= 2, where x is not equal to -1, 0 or 1 */
  unsigned i, loop_err, u_correction;
  mpfr_t p1, p2, pn, first_term, second_term;
  mpfr_prec_t res_prec, realprec, x_prec, delta, test_prec, err, lost_bits;
  MPFR_GROUP_DECL (group);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", MPFR_PREC (x), mpfr_log_prec, x, rnd_mode),
     ("legendre[%Pd]=%.*Rg", MPFR_PREC (res), mpfr_log_prec, res));

  /* first, check if x belongs to the domain [-1,1] and 0 <= n <= 2^13 = 8192.
     If it's not, res is set to NAN, and 0 is returned */
  is_within_domain &= mpfr_lessequal_p(x, __gmpfr_one);
  is_within_domain &= mpfr_greaterequal_p (x, __gmpfr_mone);

  if (!is_within_domain || n < 0 || n > 8192)
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

  res_prec = MPFR_PREC(res);
  x_prec = MPFR_PREC(x);
  /*  */
  loop_err = 14 * (n - 1);
  delta = MPFR_INT_CEIL_LOG2(res_prec) + loop_err;
  realprec = res_prec + delta;
  if (x_prec > realprec) {
    realprec = x_prec;
  }
  /* err = ceil(log2(loop_err)) + small margin (5 bits). This value is
     incremented in case of severe cancellation */
  err = MPFR_INT_CEIL_LOG2(loop_err) + 5;

  MPFR_GROUP_INIT_5 (group, realprec,
                     p1, p2, pn, first_term, second_term);

  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
_start:
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

          /* we check if a severe cancellation may occur *before* performing
             the subtraction. If the number of canceled bits is greater than
             the working precision minus the result precision, we need to
             increase the working precision */
          lost_bits = likely_cancellation (first_term, second_term, realprec);
          if (lost_bits > (realprec - res_prec))
            {
              /* increase precision by the number of lost bits
                 (due to cancellation), plus a small safety margin */
              realprec += lost_bits + 2;
              /* the error estimate is updated as well */
              err += lost_bits;

              MPFR_GROUP_REPREC_5 (group, realprec,
                                   p1, p2, pn, first_term, second_term);

              /* in case of severe cancellation, a new precision is set
                 and we go back to the beginning of the loop, so that each
                 term is recomputed with the new precision */
              goto _start;
            }

          /* pn = (first_term - second_term) / i */
          mpfr_sub (pn, first_term, second_term, MPFR_RNDN);
          mpfr_div_ui (pn, pn, i, MPFR_RNDN);

          /* p2 = p1, p1 = pn */
          mpfr_set (p2, p1, MPFR_RNDN);
          mpfr_set (p1, pn, MPFR_RNDN);
          i++;
        }

      /* we have (eventually) added the lost bits due to cancellation
         above, so test_prec is the effective precision we can use */
      test_prec = realprec - err;

      // if (mpfr_min_prec (pn) < test_prec - 1)
      //   break;

      if (MPFR_LIKELY (MPFR_CAN_ROUND (pn, test_prec, res_prec, rnd_mode)))
        break;

      if (mpfr_equal_p(pn, p1))
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
