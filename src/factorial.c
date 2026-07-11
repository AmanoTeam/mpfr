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
static const unsigned long mpfr_fac_table[] = {
    1UL                      /*  0! */
  , 1UL                      /*  1! */
  , 2UL                      /*  2! */
  , 6UL                      /*  3! */
  , 24UL                     /*  4! */
  , 120UL                    /*  5! */
  , 720UL                    /*  6! */
  , 5040UL                   /*  7! */
  , 40320UL                  /*  8! */
  , 362880UL                 /*  9! */
  , 3628800UL                /* 10! */
  , 39916800UL               /* 11! */
  , 479001600UL              /* 12! */
#if ULONG_MAX >= 2432902008176640000UL
  , 6227020800UL             /* 13! */
  , 87178291200UL            /* 14! */
  , 1307674368000UL          /* 15! */
  , 20922789888000UL         /* 16! */
  , 355687428096000UL        /* 17! */
  , 6402373705728000UL       /* 18! */
  , 121645100408832000UL     /* 19! */
  , 2432902008176640000UL    /* 20! */
#endif
};

/* when x is too large for x! to fit in an unsigned long, the product
   (numberof mpfr_fac_table) * ... * x is computed by grouping several
   consecutive integers and multiplying them together as native unsigned
   long integers, so that each group needs a single mpfr_mul_ui.
   mpfr_fac_group[b] is the maximum number of integers of b bits, i.e. in
   [2^(b-1), 2^b - 1], whose product is guaranteed to fit in an unsigned
   long: it is floor(B / b) where B is the number of bits of an unsigned
   long, since the product of g such integers is less than 2^(g*b) <= 2^B.
   The table is hardcoded for a 32-bit and a 64-bit unsigned long; entries
   beyond the table (larger b) are taken to be 1. For an unsigned long wider
   than 64 bits, the 64-bit table is used, which is still safe since the
   product is then less than 2^64 <= ULONG_MAX. */
#if ULONG_MAX >= 2432902008176640000UL  /* unsigned long has >= 64 bits */
static const unsigned char mpfr_fac_group[] = {
  0, 64, 32, 21, 16, 12, 10,  9,  8, /* b =  0.. 8 */
  7,  6,  5,  5,  4,  4,  4,  4,     /* b =  9..16 */
  3,  3,  3,  3,  3,  2,  2,  2,     /* b = 17..24 */
  2,  2,  2,  2,  2,  2,  2,  2      /* b = 25..32 */
};
#else                                   /* unsigned long has >= 32 bits */
static const unsigned char mpfr_fac_group[] = {
  0, 32, 16, 10,  8,  6,  5,  4,  4, /* b =  0.. 8 */
  3,  3,  2,  2,  2,  2,  2,  2      /* b =  9..16 */
};
#endif

/* Returns a lower bound of floor(log2(n!)).
   The exact value is log2(n!) = lgamma(n+1)/log(2) since n! = Gamma(n+1).
   A lower bound is obtained by rounding lgamma(n+1) down (MPFR_RNDD) and
   log(2) up (MPFR_RNDU), so that their quotient is rounded down too.
   The main purpose of this function is to detect overflows without computing
   the factorial.  For more "fine-grained" overflows, i.e. those in proximity
   of EXP_MAX, we handle those in the Ziv main loop.
   Note: for very large n (above MPFR_FAC_OVERFLOW_N), overflow is detected
   via the hard-coded threshold in mpfr_fac_ui, so this function is only
   called for n up to that threshold. */
static mpfr_exp_t
magnitude (unsigned long n)
{
  int inex;
  mpfr_exp_t ret;
  mpfr_prec_t realprec;
  mpfr_t lb, fn, ln2;
  MPFR_ZIV_DECL (loop);
  MPFR_GROUP_DECL (group);

  realprec = 64;

  MPFR_GROUP_INIT_3 (group, realprec, lb, fn, ln2);

  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
      /* fn = n+1, exact for n < 2^64 */
      inex = mpfr_set_ui (fn, n + 1, MPFR_RNDN);
      /* lb = lgamma(n+1) rounded down */
      inex |= mpfr_lngamma (lb, fn, MPFR_RNDD);
      /* ln2 = log(2) rounded up */
      inex |= mpfr_log (ln2, __gmpfr_two, MPFR_RNDU);
      /* lb = log2(n!) rounded down */
      inex |= mpfr_div (lb, lb, ln2, MPFR_RNDD);
      inex |= mpfr_floor (lb, lb);

      if (inex == 0)
        break;
      MPFR_ZIV_NEXT (loop, realprec);
      MPFR_GROUP_REPREC_3 (group, realprec, lb, fn, ln2);
    }

  MPFR_ZIV_FREE (loop);
  ret = mpfr_get_si (lb, MPFR_RNDN);

  MPFR_GROUP_CLEAR (group);

  return ret;
}

/* MPFR_FAC_OVERFLOW_N: smallest n such that the Stirling lower bound
     log2(n!) >= (n*log(n) - n + 1) / log(2)
   exceeds MPFR_EMAX_MAX, the maximum allowed exponent. For n greater
   or equalat this threshold, n! always overflows.
   The values are:
   - 64-bit unsigned long: MPFR_EMAX_MAX = 2^62-1; true threshold ~ 8.42e16;
     we use 9e16 = 90000000000000000 for safety;
   - 32-bit unsigned long: MPFR_EMAX_MAX = 2^30-1; true threshold ~ 4.5e7;
     we use 5e7 = 50000000 for safety */
#if ULONG_MAX > 0xFFFFFFFFUL
# define MPFR_FAC_OVERFLOW_N 90000000000000000UL  /* 9e16, for 64-bit */
#else
# define MPFR_FAC_OVERFLOW_N 50000000UL           /* 5e7,  for 32-bit */
#endif

static int
factorial (mpfr_t t, unsigned long int x, mpfr_rnd_t rnd)
{
  int inexact = 0;
  unsigned long int i;

  /* multiply t by start * (start+1) * ... * x, where start is
     numberof_const (mpfr_fac_table). Consecutive integers are grouped
     and multiplied together as native unsigned long integers, so each
     group needs a single mpfr_mul_ui */
  for (i = numberof_const (mpfr_fac_table); i <= x ; )
    {
      unsigned long p, cnt, k;
      int b, round;

      /* number of bits of i: b such that 2^(b-1) <= i < 2^b */
      b = MPFR_INT_CEIL_LOG2 (i + 1);

      /* maximum number of b-bit integers whose product fits
         in an unsigned long (1 beyond the hardcoded table) */
      if (b < (int) numberof_const (mpfr_fac_group))
        cnt = mpfr_fac_group[b];
      else
        {
          /* b is beyond the grouping table (i >= 2^16 for 32-bit, i >= 2^32
             for 64-bit): even a product of two such integers could reach
             2^(2b) > ULONG_MAX, so no grouping is safe and each integer
             must be multiplied individually. */
          cnt = 1;
        }

      /* keep all grouped integers on b bits, i.e. do not
         cross the 2^b boundary, so that the product of cnt
         of them is guaranteed to fit in an unsigned long. */
      if (b < (int) (sizeof (unsigned long) * CHAR_BIT))
        {
          unsigned long room = ((unsigned long) 1 << b) - i;
          if (cnt > room)
            cnt = room;
        }

      /* do not go past x. */
      if (cnt > x - i + 1)
        cnt = x - i + 1;

      /* p = i * (i+1) * ... * (i+cnt-1), on native integers */
      p = i;
      for (k = 1; k < cnt; k++)
        p *= i + k;
      i += cnt;

      round = mpfr_mul_ui (t, t, p, rnd);

      /* assume the first inexact product gives the sign
         of difference: is that always correct? */
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
  unsigned long int start_prod;
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
  if (MPFR_UNLIKELY (x < numberof_const (mpfr_fac_table)))
    return mpfr_set_ui (y, mpfr_fac_table[x], rnd_mode);
  else
    start_prod = mpfr_fac_table[numberof_const (mpfr_fac_table) - 1];

  /* For very large x, x! overflows for any valid emax (including the
     maximum MPFR_EMAX_MAX) */
  if (MPFR_UNLIKELY (x >= MPFR_FAC_OVERFLOW_N))
    return mpfr_overflow (y, rnd_mode, 1);

  MPFR_SAVE_EXPO_MARK (expo);

  /* Once x is at or beyond the grouping table limit, every remaining
     integer requires its own mpfr_mul_ui call, since there's no
     possible grouping. We check for overflow using
       log2(n!) = floor(lgamma(n+1)/log(2)) <= lgamma(n+1)/log(2).
     If log2(n!) cannot fit in emax, it's going to be an overflow */
  if (x >= ((unsigned long) 1 << (numberof_const (mpfr_fac_group) - 1))
      && magnitude (x) > emax)
    {
      MPFR_SAVE_EXPO_FREE (expo);
      return mpfr_overflow (y, rnd_mode, 1);
    }

  /* Initialisation of the Precision */
  Ny = MPFR_PREC (y);

  /* compute the size of intermediary variable */
  Nt = Ny + 2 * MPFR_INT_CEIL_LOG2 (x) + 7;

  mpfr_init2 (t, Nt); /* initialize of intermediary variable */

  rnd = MPFR_RNDZ;
  MPFR_ZIV_INIT (loop, Nt);
  for (;;)
    {
      MPFR_BLOCK_DECL (flags);

      /* Start the product from the largest factorial that fits in an
         unsigned long, then multiply the remaining integers below. */
      inexact = mpfr_set_ui (t, start_prod, rnd);

      MPFR_BLOCK (flags, inexact = factorial (t, x, rnd_mode));

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
