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
static mpfr_eexp_t
magnitude (unsigned long n)
{
  mpfr_eexp_t ret;
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

  ret = mpfr_get_exp_t (lb, MPFR_RNDD);

  mpfr_clear (lb);
  mpfr_clear (fn);
  mpfr_clear (ln2);

  return ret;
}

static int
factorial (mpfr_t t, unsigned long int n)
{
  int inexact;
  unsigned long int i, imax;
  unsigned int b;

  i = numberof_const (mpfr_fac_group) - 1;
  MPFR_ASSERTD (n > i);  /* n <= i handled in mpfr_fac_ui() */
  inexact = mpfr_set_ui (t, mpfr_fac_group[i], MPFR_RNDZ);

  b = MPFR_INT_CEIL_LOG2 (i + 1);
  imax = ((unsigned long) 1 << b) - 1;

  /* Multiply t by (i+1) * (i+2) * ... * n. */
  for (;;)
    {
      unsigned long cnt;

      /* Maximum number of b-bit integers whose product fits in an
         unsigned long; this is 1 once b > ULSIZE / 2, i.e. no grouping
         is possible and each integer is multiplied individually. */
      cnt = ULSIZE / b;

      if (imax > n)
        imax = n;

      while (i < imax)
        {
          unsigned long imax2, p;

          imax2 = MIN (i + cnt, imax);

          /* p = (i+1) * (i+2) * ... * imax2, on native integers */
          p = ++i;
          while (i < imax2)
            p *= ++i;

          inexact |= mpfr_mul_ui (t, t, p, MPFR_RNDZ);

          /* an overflow of an intermediate product is a real overflow: it
             occurs in the maximal exponent range (set by
             MPFR_SAVE_EXPO_MARK) and does not depend on the working
             precision Nt, so we can stop as soon as we detect one */
          if (MPFR_UNLIKELY (MPFR_BLOCK_EXCEP))
            return inexact;
        }

      if (i == n)
        return inexact;

      b++;
      imax = (imax << 1) + 1;
    }

  return inexact;
}

int
mpfr_fac_ui (mpfr_ptr y, unsigned long int n, mpfr_rnd_t rnd_mode)
{
  int inexact;
  mpfr_exp_t emax;
  mpfr_t t;         /* Variable of Intermediary Calculation */
  mpfr_prec_t Ny;   /* Precision of output variable */
  mpfr_prec_t Nt;   /* Precision of Intermediary Calculation variable */
  mpfr_prec_t err;  /* Precision of error */
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  emax = mpfr_get_emax ();

  /* For n such that n! fits in an unsigned long, we directly set y from
     the precomputed value, avoiding the costly loop of mpfr_mul_ui calls. */
  if (MPFR_UNLIKELY (n < numberof_const (mpfr_fac_group)))
    return mpfr_set_ui (y, mpfr_fac_group[n], rnd_mode);

  /* For very large n, n! overflows for any valid emax (including the
     maximum MPFR_EMAX_MAX) */
  if (MPFR_FAC_OVERFLOW_N > 0 && MPFR_UNLIKELY (n >= MPFR_FAC_OVERFLOW_N))
    return mpfr_overflow (y, rnd_mode, 1);

  MPFR_SAVE_EXPO_MARK (expo);

  /* Once n >= MPFR_FAC_NO_GROUPING_N, every remaining integer requires
     its own multiplication, so it is worth checking for overflow using
       log2(n!) = floor(lgamma(n+1)/log(2)) <= lgamma(n+1)/log(2).
     If log2(n!) > emax, it's going to be an overflow. */
  if (n >= MPFR_FAC_NO_GROUPING_N && magnitude (n) > emax)
    {
      MPFR_SAVE_EXPO_FREE (expo);
      return mpfr_overflow (y, rnd_mode, 1);
    }

  Ny = MPFR_PREC (y);
  Nt = Ny + 2 * MPFR_INT_CEIL_LOG2 (n) + 7;

  mpfr_init2 (t, Nt);

  MPFR_ZIV_INIT (loop, Nt);
  for (;;)
    {
      MPFR_BLOCK_DECL (flags);

      MPFR_BLOCK (flags, inexact = factorial (t, n));

      /* Since we rounded toward zero (MPFR_RNDZ), an intermediate overflow
         necessarily is a real overflow. And once the working precision is
         large enough, an overflow will necessarily be detected. */
      if (MPFR_UNLIKELY (MPFR_OVERFLOW (flags)))
        {
          MPFR_ZIV_FREE (loop);
          mpfr_clear (t);
          MPFR_SAVE_EXPO_FREE (expo);
          return mpfr_overflow (y, rnd_mode, 1);
        }

      /* We have an error bound expressed with a factor of the typical
         form ku/(1-ku), where k is the number of inexact products and
         u = 2^(1-p). And thanks to a "MPFR_INT_CEIL_LOG2 (n)" term in
         the initial working precision (and even more), p > 2+log2(k),
         so that ku < 1/2. So the factor is less than 2ku, and the number
         of lost bits is less than log2(2k) < 1 + MPFR_INT_CEIL_LOG2 (n).
         Note: One could do better, but this should not be noticeable in
         practice, because it is expected that even with this bound, the
         MPFR_CAN_ROUND will succeed in general. */
      err = Nt - 1 - MPFR_INT_CEIL_LOG2 (n);

      if (MPFR_LIKELY (!inexact || MPFR_CAN_ROUND (t, err, Ny, rnd_mode)))
        break;

      MPFR_ZIV_NEXT (loop, Nt);
      mpfr_set_prec (t, Nt);
    }
  MPFR_ZIV_FREE (loop);

  inexact = mpfr_set (y, t, rnd_mode);

  mpfr_clear (t);
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inexact, rnd_mode);
}
