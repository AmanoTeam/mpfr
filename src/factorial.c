/* mpfr_fac_ui -- factorial of a non-negative integer

Copyright 2001, 2004-2026 Free Software Foundation, Inc.
Contributed by the Pascaline and Caramba projects, INRIA.

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

/* table of the values of n! that fit in an unsigned long. For a 32-bit
   unsigned long, n! fits up to 12! = 479001600, and for a 64-bit unsigned
   long, up to 20! = 2432902008176640000 */
static const unsigned long mpfr_fac_group[] = {
    1UL                   /*  0! */
  , 1UL
  , 2UL
  , 6UL
  , 24UL
  , 120UL
  , 720UL
  , 5040UL
  , 40320UL
  , 362880UL
  , 3628800UL
  , 39916800UL
  , 479001600UL
#if ULONG_MAX >> 31 >> 31 != 0
  , 6227020800UL
  , 87178291200UL
  , 1307674368000UL
  , 20922789888000UL
  , 355687428096000UL
  , 6402373705728000UL
  , 121645100408832000UL
  , 2432902008176640000UL /* 20! */
#endif /* ULONG_MAX >> 31 >> 31 != 0 */
};

/* MPFR_FAC_OVERFLOW_N is a threshold on n above which n! is guaranteed to
   overflow for any valid exponent range, so that mpfr_fac_ui can return an
   overflow directly. See algorithms.tex for further details.
   MPFR_EMAX_MAX is currently the integer part of MPFR_EXP_MAX / 2. Since
     1073741825 >= 2^30 - 1
     4611686018427387924 >= 2^62 - 1
   the first condition below is satisfied with a 32-bit exponent and the
   second condition below is satisfied with a 64-bit exponent.
   But the code below remains valid even if the MPFR_EMAX_MAX definition
   is changed (MPFR_EMAX_MAX is compared with values obtained thanks to a
   weakened form of Stirling's bound formula, see algorithms.tex).
*/
#if MPFR_EMAX_MAX <= 1073741825
# define MPFR_FAC_OVERFLOW_N 44787928UL
#elif ULONG_MAX >> 31 >> 31 != 0 && MPFR_EMAX_MAX <= 4611686018427387924
# define MPFR_FAC_OVERFLOW_N 84182992257887725UL
#else
/* Probably a 128-bit exponent or a non-standard configuration. */
# define MPFR_FAC_OVERFLOW_N 0 /* no early overflows detection */
#endif

/* number of bits of an unsigned long, assuming no padding bits */
#define ULSIZE (sizeof (unsigned long) * CHAR_BIT)

/* threshold on n above which no grouping is possible: the number of b-bit
   integers whose product fits in an unsigned long, floor (ULSIZE / b), is 1
   as soon as b > ULSIZE / 2, so each integer >= 2^(ULSIZE/2) needs its own
   mpfr_mul_ui call */
#define MPFR_FAC_NO_GROUPING_N ((unsigned long) 1 << (ULSIZE / 2))

/* the main purpose of this function is to detect overflows without computing
   the factorial, by returning a lower bound of floor(log2(n!)). The exact
   value is log2(n!) = lgamma(n+1)/log(2), and a lower bound is obtained by
   rounding lgamma(n+1) down (MPFR_RNDD) and log(2) up (MPFR_RNDU), so that
   their quotient is rounded down too.
   Note: for very large n (above MPFR_FAC_OVERFLOW_N), overflow is detected
   via the hard-coded threshold in mpfr_fac_ui, so this function is only
   called for n up to that threshold */
static mpfr_exp_t
magnitude (unsigned long n)
{
  mpfr_exp_t ret;
  mpfr_t lb, fn, ln2;

  /* we check that n < ULONG_MAX, so n+1 does not overflow */
  MPFR_ASSERTD (n < ULONG_MAX);

  mpfr_init2 (lb, 64);
  mpfr_init2 (fn, 64);
  mpfr_init2 (ln2, 64);

  mpfr_set_ui (fn, n + 1, MPFR_RNDD);
  mpfr_lngamma (lb, fn, MPFR_RNDD);
  mpfr_log (ln2, __gmpfr_two, MPFR_RNDU);
  mpfr_div (lb, lb, ln2, MPFR_RNDD);
  mpfr_floor (lb, lb);

  ret = mpfr_get_si (lb, MPFR_RNDD);

  mpfr_clear (lb);
  mpfr_clear (fn);
  mpfr_clear (ln2);

  return ret;
}

static int
factorial (mpfr_t t, unsigned long int x, mpfr_rnd_t rnd)
{
  int inexact = 0;
  unsigned long int i;

  /* multiply t by 2 * 3 * ... * x. Consecutive integers are grouped
     and multiplied together as native unsigned long integers, so each
     group needs a single mpfr_mul_ui. The maximum number of integers
     of b bits, i.e. in [2^(b-1), 2^b - 1], whose product is guaranteed
     to fit in an unsigned long, is floor(ULSIZE / b), since the product
     of g such integers is less than 2^(g*b) <= 2^ULSIZE */
  for (i = 2; i <= x ; )
    {
      unsigned long p, cnt, room, j;
      int b, round;

      b = MPFR_INT_CEIL_LOG2 (i + 1);

      MPFR_ASSERTD (b >= 1 && b < ULSIZE);

      /* maximum number of b-bit integers whose product fits in an unsigned
         long; this is 1 once b > ULSIZE / 2, i.e. no grouping is possible
         and each integer is multiplied individually */
      cnt = ULSIZE / b;

      /* keep all grouped integers on b bits, i.e. do not
         cross the 2^b boundary, so that the product of cnt
         of them is guaranteed to fit in an unsigned long.
         Note that room >= 1 since i < 2^b, thus cnt >= 1 */
      room = ((unsigned long) 1 << b) - i;
      if (cnt > room)
        cnt = room;

      /* p = i * (i+1) * ... * min (i+cnt-1, x), on native integers;
         the j <= x condition also stops the group at x */
      p = 1;
      for (j = i; j < i + cnt && j <= x; j++)
        p *= j;
      i = j;

      round = mpfr_mul_ui (t, t, p, rnd);

      /* assume the first inexact product gives the sign
         of difference: is that always correct?
         FIXME: no. With a precision of 4 bits, if we approximate 7!
         with rounding to nearest, we get successively 1, 2, 6, 24,
         120, then 120*6 is rounded to 704, and 704*7 is rounded to 5120.
         The first inexact product is 120*6 which is smaller than 720,
         but the final result 5120 is larger than 7!=5040. */
      if (inexact == 0)
        inexact = round;

      /* an overflow of an intermediate product is a real overflow: it
         occurs in the maximal exponent range (set by
         MPFR_SAVE_EXPO_MARK) and does not depend on the working
         precision Nt, so we can stop as soon as we detect one */
      if (MPFR_UNLIKELY (MPFR_BLOCK_EXCEP))
        break;
    }

  return inexact;
}

int
mpfr_fac_ui (mpfr_ptr y, unsigned long int x, mpfr_rnd_t rnd_mode)
{
  int inexact;
  mpfr_exp_t emax;
  mpfr_t t;         /* Variable of Intermediary Calculation */
  mpfr_prec_t Ny;   /* Precision of output variable */
  mpfr_prec_t Nt;   /* Precision of Intermediary Calculation variable */
  mpfr_prec_t err;  /* Precision of error */
  mpfr_rnd_t rnd;
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  emax = mpfr_get_emax ();

  /* for x such that x! fits in an unsigned long, we directly set y from
     the hardcoded value, avoiding the costly loop of mpfr_mul_ui calls */
  if (MPFR_UNLIKELY (x < numberof_const (mpfr_fac_group)))
    return mpfr_set_ui (y, mpfr_fac_group[x], rnd_mode);

  /* for very large x, x! overflows for any valid emax (including the
     maximum MPFR_EMAX_MAX) */
  if (MPFR_FAC_OVERFLOW_N > 0 && MPFR_UNLIKELY (x >= MPFR_FAC_OVERFLOW_N))
    return mpfr_overflow (y, rnd_mode, 1);

  MPFR_SAVE_EXPO_MARK (expo);

  /* once x >= MPFR_FAC_NO_GROUPING_N, every remaining integer requires its
     own mpfr_mul_ui call, so it is worth checking for overflow using
       log2(n!) = floor(lgamma(n+1)/log(2)) <= lgamma(n+1)/log(2).
     If log2(n!) cannot fit in emax, it's going to be an overflow */
  if (x >= MPFR_FAC_NO_GROUPING_N
      && magnitude (x) > emax)
    {
      MPFR_SAVE_EXPO_FREE (expo);
      return mpfr_overflow (y, rnd_mode, 1);
    }

  Ny = MPFR_PREC (y);
  Nt = Ny + 2 * MPFR_INT_CEIL_LOG2 (x) + 7;

  mpfr_init2 (t, Nt);

  rnd = MPFR_RNDZ;
  MPFR_ZIV_INIT (loop, Nt);
  for (;;)
    {
      MPFR_BLOCK_DECL (flags);

      /* t = 1 is exact whatever the precision, so that the ternary value
         of the whole product is the one returned by factorial() */
      mpfr_set_ui (t, 1, rnd);

      MPFR_BLOCK (flags, inexact = factorial (t, x, rnd));

      /* since x! > 0, the overflow always yields +Inf */
      if (MPFR_UNLIKELY (MPFR_OVERFLOW (flags)))
        {
          MPFR_ZIV_FREE (loop);
          mpfr_clear (t);
          MPFR_SAVE_EXPO_FREE (expo);
          return mpfr_overflow (y, rnd_mode, 1);
        }

      err = Nt - 1 - MPFR_INT_CEIL_LOG2 (Nt);

      if (MPFR_LIKELY (!inexact || MPFR_CAN_ROUND (t, err, Ny, rnd_mode)))
        {
          /* If inexact = 0, then t is exactly x!, so round is the
             correct inexact flag.
             Otherwise, t != x! since we rounded to zero or away. */
          int round = mpfr_set (y, t, rnd_mode);
          if (inexact == 0)
            {
              inexact = round;
              break;
            }
          else if ((inexact < 0 && round <= 0) ||
                   (inexact > 0 && round >= 0))
            break;
          else /* inexact and round have opposite signs: we cannot
                  compute the inexact flag. Restart using the
                  symmetric rounding. */
            rnd = (rnd == MPFR_RNDZ) ? MPFR_RNDU : MPFR_RNDZ;
        }
      MPFR_ZIV_NEXT (loop, Nt);
      mpfr_set_prec (t, Nt);
    }
  MPFR_ZIV_FREE (loop);

  mpfr_clear (t);
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inexact, rnd_mode);
}
