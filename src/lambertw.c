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

#define W0_N_POLY      16
#define W0_POLY_DEGREE  8

/* in this table we store the coefficients of the 8-degree polynomials
   used as the initial guess to compute W0. See algorithms.tex for more
   information on how they are calculated */
static const double w0_poly_coeffs_range[W0_N_POLY][W0_POLY_DEGREE+1] = {
  { -9.544845049527667e-5, 0.00018474505346700832, -0.0001283213366529374,
    0.0002934622830968015, -0.0009249200421781934, 0.002390169033275487,
    -0.007618699877818394, 0.045176130353980194, -0.8453614329874922 },
  { -2.3659052141023946e-6, 7.4058565308182645e-6, -1.798442855079289e-5,
    6.19306131522831e-5, -0.0002334951745950057, 0.0009584790471684754,
    -0.004870812656218898, 0.04369493424052037, -0.7562369080488668 },
  { -0.0005478745121648273, 0.0009951850887610979, -0.000482623907995086,
    0.0010771622292033664, -0.0037421302971768367, 0.009001890002177297,
    -0.026066643075838555, 0.1274723453753279, -0.5472529619674287 },
  { -0.00040087402340246027, 0.0007998794499253494, -0.0006342796473689231,
    0.0014646955003262104, -0.004468884999893969, 0.011613382585588142,
    -0.03594225267984987, 0.17490687647478628, -0.2591711018190737 },
  { -0.0004697249710214173, 0.0009777285457924596, -0.0009097687952145921,
    0.0021431711067555316, -0.0063637567420890395, 0.016832365747407464,
    -0.05200482786876532, 0.23847301437678164, 0.20634179899486457 },
  { -0.0006920326810875182, 0.0014439456077716762, -0.001348453309862649,
    0.003147027489058536, -0.009232901123789574, 0.023917120251567277,
    -0.07116774736316217, 0.30305970999878734, 0.8190289382195982 },
  { -0.0009909398347625856, 0.002049427136446889, -0.0018451550881029076,
    0.004236922973801338, -0.012348836913655122, 0.03116775739052421,
    -0.0893265693137629, 0.35859757535450737, 1.5687187910944447 },
  { -0.00129394911422332, 0.002650161083060788, -0.0022931136753530377,
    0.005192987294889481, -0.015102969562114273, 0.037345079824128914,
    -0.10409799188305148, 0.40149115277417186, 2.429726378010434 },
  { -0.001549180298107751, 0.0031496279220049267, -0.002643412936250887,
    0.005929368478673035, -0.017244895897477325, 0.042059619163970476,
    -0.11512493152296555, 0.4329021539159655, 3.3742209248873585 },
  { -0.001742326701356674, 0.0035252537413985267, -0.0028990911523159813,
    0.006463779618263067, -0.01881000844308342, 0.045483284701620906,
    -0.12308839737122172, 0.4556137128974086, 4.37932512833984 },
  { -0.001881989233314381, 0.0037964199542476233, -0.0030823165743224744,
    0.006846783715612749, -0.019935147822805803, 0.04794774731641583,
    -0.12884371749486456, 0.47222256644307403, 5.4284280241852105 },
  { -0.0019824130149203555, 0.00399158532353859, -0.0032149234524559527,
    0.007124804652124798, -0.020752263929593055, 0.04974643357358327,
    -0.13307957896755107, 0.48464287893499824, 6.510033252338019 },
  { -0.0020557546791875442, 0.004134421567372127, -0.003313071561030985,
    0.007331412922698941, -0.02135901592622225, 0.05109038830635446,
    -0.1362752526525268, 0.49416879477018144, 7.616237887843317 },
  { -0.002110627609332965, 0.004241552908693455, -0.0033876249348644117,
    0.0074890024691269905, -0.021821266459289737, 0.0521206692117194,
    -0.13874812080496352, 0.5016542282762467, 8.741532903265194 },
  { -0.0021527542077419766, 0.004323996237115148, -0.003445696377254052,
    0.007612220966392506, -0.022182277769603437, 0.05292990772470599,
    -0.14070695780211315, 0.5076654781053735, 9.88199000551283 },
  { -0.0021858862926890034, 0.004388977490341032, -0.0034919644021215823,
    0.007710724494479912, -0.022470590009471738, 0.05357944427792002,
    -0.1422909411721637, 0.5125850933995791, 11.034737543928468 }};

/* here we have the pre-calculated upper bounds of the interval [a,b] of each
   polynomial. See algorithms.tex for further information */
static const double w0_range_bounds[W0_N_POLY+1] = {
  -0.368, -0.36, -0.35, -0.3, -0.1, 0.60725964, 3.1083407, 11.952909, 43.22994,
  153.83485, 544.9668, 1928.126, 6819.3895, 24116.356, 85283.587, 301589.16,
  1066510.2 };

/* Horner method for polynomial evaluation */
static int
poly_horner (mpfr_t res, const mpfr_t *c, size_t ncoeff, const mpfr_t x,
             mpfr_rnd_t rnd)
{
  int t;
  size_t i;

  if (ncoeff == 0)
    {
      mpfr_set_zero (res, +1);
      return 0;
    }

  t = mpfr_set (res, c[0], rnd);

  for (i = 1; i < ncoeff; i++)
    t |= mpfr_fma (res, res, x, c[i], rnd);

  return t;
}

static int
eval_poly (mpfr_ptr res, const double *coeffs, size_t n_coeffs,
           mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int ternary;
  size_t i;
  mpfr_t *c;

  c = (mpfr_t *) malloc (n_coeffs * sizeof (mpfr_t));
  for (i = 0; i < n_coeffs; i++)
    {
      mpfr_init2 (c[i], MPFR_PREC (res));
      mpfr_set_d (c[i], coeffs[i], MPFR_RNDN);
    }

  ternary = poly_horner (res, c, n_coeffs, x, rnd_mode);

  for (i = 0; i < n_coeffs; i++)
    {
      mpfr_clear (c[i]);
    }
  free (c);

  return ternary;
}

#define EVAL_POLY_RANGE(res,k,x,rnd) \
  eval_poly (res, w0_poly_coeffs_range[k], W0_POLY_DEGREE + 1, x, rnd)

/* the correctly rounded implementation of the constant 1/e, where e is
   Euler'n number \approx 2.7182 */
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
      /* W_0(-1/e) = W_{-1}(-1/e) = -1 */
      return mpfr_set (res, __gmpfr_mone, rnd_mode);
    }

  /* the leftmost point of the real domain of both W_0 and w_{-1} is x = -1/e,
     so we set res = NaN */
  MPFR_SET_NAN (res);
  MPFR_RET_NAN;
}

/* we have to map back the domain [a,b] to [-1,1], which is the domain of the
   interpolating polynomials, by (2*x-a-b)/(b-a). For more details concerning
   the interpolating polynomials and their domain, see algorithms.tex */
static void
map_x_in_ab (mpfr_t t, mpfr_srcptr x, double a, double b)
{
  mpfr_mul_2ui (t, x, 1, MPFR_RNDN);
  mpfr_sub_d (t, t, a + b, MPFR_RNDN);
  mpfr_div_d (t, t, b - a, MPFR_RNDN);
}

/* sets res to an initial approximation of W0(x), used as the starting point
   of the Halley iteration. The domain [-1/e,+Inf) is split into subregions,
   on each of which W0 is approximated by a polynomial; see algorithms.tex for
   how the subregions and the coefficients are computed */
static void
initial_guess_w0 (mpfr_ptr res, mpfr_srcptr x, mpfr_srcptr inve)
{
  int k = 0, is_large;
  double a, b, n;
  mpfr_t t, m;

  mpfr_init2 (m, 53);
  mpfr_set_d (m, w0_range_bounds[W0_N_POLY], MPFR_RNDN);

  is_large = mpfr_cmp (x, m) >= 0;

  mpfr_clear (m);

  if (is_large)
    {
      /* TODO: handle large x using the expansion:
         W_0(x) \approx L_1 - L_2 + L_2/L_1 with
           L_1 = log x;
           L_2 = log L_1.
         for now we return 11.442562431228513, that is W0(1066510.2) */
      mpfr_set_d(res, 11.442562431228513, MPFR_RNDN);
      return;
    }

  /* mpfr_get_d is not going to overflow, since at this point -1/e < n
     < 1066510.2, so we can store x in a double */
  n = mpfr_get_d (x, MPFR_RNDN);
  if (n > w0_range_bounds[0])
    {
      for (; k < W0_N_POLY && n >= w0_range_bounds[k + 1]; k++);
      a = w0_range_bounds[k];
      b = w0_range_bounds[k+1];
    }
  else
    {
      k = 0;
      a = w0_range_bounds[0];
      b = w0_range_bounds[1];
    }

  mpfr_init2 (t, MPFR_PREC (res));

  map_x_in_ab (t, x, a, b);
  EVAL_POLY_RANGE (res, k, t, MPFR_RNDN);

  mpfr_clear (t);
}

int
mpfr_lambertwm1 (mpfr_ptr res, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int cmp_boundary;
  mpfr_t inve;

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", MPFR_PREC (x), mpfr_log_prec, x, rnd_mode),
     ("lambertwm1[%Pd]=%.*Rg", MPFR_PREC (res), mpfr_log_prec, res));

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
  int cmp_boundary, ternary;
  mpfr_t inve, w, ew, f, dw, denom, t;
  mpfr_prec_t realprec;
  mpfr_exp_t prev_dw_exp = MPFR_EXP_MAX;
  MPFR_GROUP_DECL (group);

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

  if (cmp_boundary <= 0)
    {
      mpfr_clear (inve);
      return early_exit_on_boundary (res, cmp_boundary, rnd_mode);
    }

  realprec = MPFR_PREC (res) + 20;

  MPFR_GROUP_INIT_6 (group, realprec, w, ew, f, dw, denom, t);

  initial_guess_w0 (w, x, inve);

  mpfr_clear (inve);

  for (;;)
    {
      mpfr_exp_t dw_exp, d;

      mpfr_exp (ew, w, MPFR_RNDN);
      mpfr_fms (f, w, ew, x, MPFR_RNDN);

      /* the initial guess is refined with Halley's iteration. In the code
         below, an Halley step calculates f(w) = w*exp(w) - x. Its derivatives
          are f'(w) = (w+1)*exp(w) and f''(w) = (w+2)*exp(w). For further
          details, see algorithms.tex */
      mpfr_add_ui (t, w, 2, MPFR_RNDN);
      mpfr_mul (t, t, f, MPFR_RNDN);
      mpfr_add_ui (denom, w, 1, MPFR_RNDN);
      mpfr_div (t, t, denom, MPFR_RNDN);
      mpfr_div_2ui (t, t, 1, MPFR_RNDN);
      mpfr_mul (denom, denom, ew, MPFR_RNDN);
      mpfr_sub (denom, denom, t, MPFR_RNDN);

      if (MPFR_IS_ZERO (denom))
        {
          MPFR_LOG_MSG (("Trying to compute W0(%.*Rg); the Halley iteration "
                         "gave a zero denominator. Exit...", x));
          MPFR_SET_NAN (res);
          MPFR_RET_NAN;
        }

      mpfr_div (dw, f, denom, MPFR_RNDN);

      /* if f (and therefore dw) is zero, we already reached a fixed
         point of the Halley iteration */
      if (MPFR_IS_ZERO (dw))
        break;

      /* if dw is singular it means that something in the code above was not
         handled properly (e.g. an overflow). Nothing left to do */
      MPFR_ASSERTD (!MPFR_IS_SINGULAR (dw));

      mpfr_sub (w, w, dw, MPFR_RNDN);

      dw_exp = MPFR_GET_EXP (dw);

      /* since x != 0, w = 0 cannot be a fixed point */
      if (MPFR_IS_ZERO (w))
        {
          MPFR_LOG_MSG (("Trying to compute w = W0(%.*Rg); the Halley iteration"
                         "calculate w = 0. Exit...", x));
          MPFR_SET_NAN (res);
          MPFR_RET_NAN;
        }

      d = MPFR_GET_EXP (w) - dw_exp;

      if (d >= (mpfr_exp_t) realprec
          || (d >= (mpfr_exp_t) (realprec / 2) && dw_exp >= prev_dw_exp))
        break;
      prev_dw_exp = dw_exp;
    }

  ternary = mpfr_set (res, w, rnd_mode);

  MPFR_GROUP_CLEAR (group);

  MPFR_RET (ternary);
}
