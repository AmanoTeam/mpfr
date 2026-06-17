/* mpfr_lambertw -- computes the two real branches W_0(x) and W_{-1}(x)
   of the LambertW complex function

Copyright 2026 Free Software Foundation, Inc.
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
mpfr_const_inve (mpfr_ptr res)
{
  int inex, ternary;
  mpfr_t e, inve;
  mpfr_exp_t err;
  mpfr_prec_t realprec, res_prec;

  MPFR_ZIV_DECL (loop);
  MPFR_GROUP_DECL (group);

  res_prec = MPFR_PREC (res);
  realprec = res_prec + MPFR_INT_CEIL_LOG2 (res_prec);

  MPFR_GROUP_INIT_2 (group, realprec, e, inve);
  MPFR_ZIV_INIT (loop, realprec);

  for (;;)
    {
      /* e = exp(1), with error <= 1/2 ulp(e) since mpfr_exp is correctly
         rounded to nearest */
      inex = mpfr_exp (e, __gmpfr_one, MPFR_RNDN);
      /* inve = 1/e. By the generic error of the division (see algorithms.tex),
         with an exact numerator and a denominator e known with error
         <= 1/2 ulp(e), we get
            error(inve) <= (1/2 + 2*1*2*(1/2)) ulp(inve) = 5/2 ulp(inve)
                        <= 2^2 ulp(inve),
         hence err = 2 */
      inex |= mpfr_ui_div (inve, 1, e, MPFR_RNDN);
      err = 2;

      /* if inex = 0, the computation was exact, thus inve is exactly 1/e;
         otherwise inve approximates 1/e with error <= 2^err ulp(inve), and we
         use MPFR_CAN_ROUND to check whether inve can be rounded to res_prec */
      if (inex == 0
          || MPFR_CAN_ROUND (inve, realprec - err, res_prec, MPFR_RNDN))
        break;

      MPFR_ZIV_NEXT (loop, realprec);
      MPFR_GROUP_REPREC_2 (group, realprec, e, inve);
    }

  MPFR_ZIV_FREE (loop);

  ternary = mpfr_set (res, inve, MPFR_RNDN);

  MPFR_GROUP_CLEAR (group);

  return ternary;
}

static int
early_exit_on_boundary (mpfr_ptr res, int cmp_boundary, mpfr_rnd_t rnd_mode)
{
  if (cmp_boundary == 0)
    {
      /* W0(-1/e) = -1 */
      return mpfr_set (res, __gmpfr_mone, rnd_mode);
    }
  
  /* the leftmost point of the real domain of both W_0 and w_{-1} is x = -1/e,
     so we set res = NaN */
  MPFR_SET_NAN (res);
  MPFR_RET_NAN;
}

int
mpfr_lambertw1 (mpfr_ptr res, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int cmp_boundary;
  mpfr_t inve;

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", MPFR_PREC (x), mpfr_log_prec, x, rnd_mode),
     ("lambertw1[%Pd]=%.*Rg", MPFR_PREC (res), mpfr_log_prec, res));

  if (MPFR_IS_ZERO (x))
    {
      /* W_{-1}(0) = -Inf */
      MPFR_SET_INF (res);
      MPFR_SET_NEG (res);
      MPFR_RET (0);
    }

  if (MPFR_IS_NAN (x) || MPFR_IS_POS (x))
    {
      MPFR_SET_NAN (res);
      MPFR_RET_NAN;
    }

  /* we need -1/e correctly rounded to MPFR_PREC (x) */
  mpfr_init2 (inve, MPFR_PREC (x));
  mpfr_const_inve (inve);
  MPFR_SET_NEG (inve);
  cmp_boundary = mpfr_cmp (x, inve);
  mpfr_clear (inve);

  if (cmp_boundary <= 0)
    return early_exit_on_boundary (res, cmp_boundary, rnd_mode);

  return 0;
}

int
mpfr_lambertw0 (mpfr_ptr res, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int cmp_boundary;
  mpfr_t inve;

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", MPFR_PREC (x), mpfr_log_prec, x, rnd_mode),
     ("lambertw0[%Pd]=%.*Rg", MPFR_PREC (res), mpfr_log_prec, res));

  if (MPFR_IS_ZERO (x))
    {
      MPFR_SET_ZERO (res);
      MPFR_RET (0);
    }

  if (MPFR_IS_NAN (x))
    {
      MPFR_SET_NAN (res);
      MPFR_RET_NAN;
    }

  if (MPFR_IS_INF (x) && MPFR_IS_POS (x))
    {
      /* W_0(+Inf) = +Inf */
      MPFR_SET_INF (res);
      MPFR_SET_POS (res);
      MPFR_RET (0);
    }

  /* we need -1/e correctly rounded to MPFR_PREC (x) */
  mpfr_init2 (inve, MPFR_PREC (x));
  mpfr_const_inve (inve);
  MPFR_SET_NEG (inve);
  cmp_boundary = mpfr_cmp (x, inve);
  mpfr_clear (inve);

  if (cmp_boundary <= 0)
    return early_exit_on_boundary (res, cmp_boundary, rnd_mode);

  return 0;
}
