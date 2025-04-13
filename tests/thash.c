/* Test file for mpfr_hash

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

#include "mpfr-test.h"

#define N_SAMPLES 200
#define H_PI      126454210U
#define H_LOG2    1348199468U
#define H_EULER   3879845459U
#define H_CATALAN 887850484U

#if SIZEOF_UNSIGNED_INT == 4
#define H_FORMAT "%u"
#elif SIZEOF_UNSIGNED_LONG == 4
#define H_FORMAT "%lu"
#endif /* SIZEOF_UNSIGNED_INT */

static mpfr_t pos[N_SAMPLES], neg[N_SAMPLES];

static void
free_mpfr_arr (mpfr_t *arr, size_t size)
{
  size_t i;
  for (i = 0; i < size; i++)
    mpfr_clear (arr[i]);
}

static void
test_zero (void)
{
  int i, j;
  mpfr_t pZero, nZero;
  mpfr_digest_t pZero_hash, nZero_hash;

  mpfr_init2 (pZero, 200);
  mpfr_set_zero (pZero, 1);
  mpfr_init2 (nZero, 200);
  mpfr_set_zero (nZero, -1);

  pZero_hash = mpfr_hash (pZero);
  nZero_hash = mpfr_hash (nZero);

  /* Zeros with different signs should be equal */
  if (pZero_hash != nZero_hash)
    {
      fprintf (stderr, "hash for +0.0 and -0.0 should be equal."
                       "H(+0.0) = " H_FORMAT "; H(-0.0) = " H_FORMAT "\n",
                       pZero_hash, nZero_hash);
      exit (1);
    }

  mpfr_clear (pZero);
  mpfr_clear (nZero);

  /* Zeros with different precisions should be equal */
  for (i = 0, j = MPFR_PREC_MIN; i < N_SAMPLES; i++, j++)
    {
      mpfr_init2 (pos[i], j);
      mpfr_set_zero (pos[i], 1);
      mpfr_init2 (neg[i], j);
      mpfr_set_zero (neg[i], -1);
    }

  for (i = 0, j = 1; j < N_SAMPLES; i++, j++)
    {
      if (mpfr_hash (pos[i]) != mpfr_hash (pos[j])
        && mpfr_hash (neg[i]) != mpfr_hash (neg[j]))
      {
        fprintf (stderr, "All zeros should be hashed the same, regardless their "
                         "precision. i: %d j: %d\n", i , j);
        exit (1);
      }
    }

  free_mpfr_arr (pos, N_SAMPLES);
  free_mpfr_arr (neg, N_SAMPLES);

  mpfr_free_cache ();
}

static void
test_inf (void)
{
  int i, j;
  mpfr_t pInf, nInf;
  mpfr_digest_t pInf_hash, nInf_hash;

  mpfr_init2 (pInf, 200);
  mpfr_set_inf (pInf, 1);
  mpfr_init2 (nInf, 200);
  mpfr_set_inf (nInf, -1);

  pInf_hash = mpfr_hash (pInf);
  nInf_hash = mpfr_hash (nInf);

  /* H(+Inf) and H(-Inf) should not be equal */
  if (pInf_hash == nInf_hash)
    {
      fprintf (stderr, "hash for +Inf and -Inf should not be equal."
                       "H(+Inf) = " H_FORMAT "; H(-Inf) = " H_FORMAT "\n",
                       pInf_hash, nInf_hash);
      exit (1);
    }

  mpfr_clear (pInf);
  mpfr_clear (nInf);

  /* Infinite numbers with different precisions should be equal */
  for (i = 0, j = MPFR_PREC_MIN; i < N_SAMPLES; i++, j++)
    {
      mpfr_init2 (pos[i], j);
      mpfr_set_inf (pos[i], 1);
      mpfr_init2 (neg[i], j);
      mpfr_set_inf (neg[i], -1);
    }

  for (i = 0, j = 1; j < N_SAMPLES; i++, j++)
  {
    if (mpfr_hash (pos[i]) != mpfr_hash (pos[j]))
      {
        fprintf (stderr, "All +Inf should be hashed the same, regardless their "
                         "precision. i: %d j: %d\n", i , j);
        exit (1);
      }
    if (mpfr_hash (neg[i]) != mpfr_hash (neg[j]))
      {
        fprintf (stderr, "All -Inf should be hashed the same, regardless their "
                         "precision. i: %d j:%d\n", i , j);
        exit (1);
      }
  }

  free_mpfr_arr (pos, N_SAMPLES);
  free_mpfr_arr (neg, N_SAMPLES);

  mpfr_free_cache ();
}

static void
test_nan (void)
{
  int i, j;
  mpfr_t nan, unconventional_nan;
  mpfr_digest_t nan_hash, unconventional_nan_hash;

  mpfr_init2 (nan, 200);
  mpfr_set_nan (nan);
  mpfr_init2 (unconventional_nan, 200);
  mpfr_set_nan (unconventional_nan);

  MPFR_CHANGE_SIGN (unconventional_nan);

  nan_hash = mpfr_hash (nan);
  unconventional_nan_hash = mpfr_hash (unconventional_nan);

  /* mpfr_hash should ignore the sign of NAN */
  if (nan_hash != unconventional_nan_hash)
    {
      fprintf (stderr, "hash for +NAN and -NAN should be equal."
                       "H(+NAN) = " H_FORMAT "; H(-NAN) = " H_FORMAT "\n",
                       nan_hash, unconventional_nan_hash);
      exit (1);
    }

  mpfr_clear (nan);
  mpfr_clear (unconventional_nan);

  /* NAN numbers with different precision should be equal */
  for (i = 0, j = MPFR_PREC_MIN; i < N_SAMPLES; i++, j++)
    {
      mpfr_init2 (pos[i], j);
      mpfr_set_nan (pos[i]);
    }

  for (i = 0, j = 1; j < N_SAMPLES; i++, j++)
    {
      if (mpfr_hash (pos[i]) != mpfr_hash (pos[j]))
       {
          fprintf (stderr, "All NANs should be hashed the same, regardless their "
                           "precision. i: %d j: %d\n", i , j);
          exit (1);
       }
    }

  free_mpfr_arr (pos, N_SAMPLES);

  mpfr_free_cache ();
}

static void
test_precision (void)
{
  mpfr_t low_prec, high_prec;
  mpfr_digest_t hash_low_prec, hash_high_prec;

  mpfr_init2 (low_prec, MPFR_PREC_MIN);
  mpfr_set_d (low_prec, 1.0, MPFR_RNDD);

  /* 20 bit precision:
       - 1 limb for 32-bit archs;
       - 1 limb for 64-bit archs */
  mpfr_init2 (high_prec, 20);
  mpfr_set_d (high_prec, 1.0, MPFR_RNDD);

  hash_low_prec = mpfr_hash(low_prec);
  hash_high_prec = mpfr_hash(high_prec);

  if (hash_low_prec != hash_high_prec)
    {
      fprintf (stderr, "1.0 (prec. 1) and 1.0 (prec. 20) should have the "
                       "same hash.\n");
      exit (1);
    }

  mpfr_clear (high_prec);

  /* 50 bit precision:
       - 2 limbs for 32-bit archs;
       - 1 limb for 64-bit archs */
  mpfr_init2 (high_prec, 50);
  mpfr_set_d (high_prec, 1.0, MPFR_RNDD);

  hash_low_prec = mpfr_hash(low_prec);
  hash_high_prec = mpfr_hash(high_prec);

  if (hash_low_prec != hash_high_prec)
    {
      fprintf (stderr, "1.0 (prec. 1) and 1.0 (prec. 50) should have the "
                       "same hash.\n");
      exit (1);
    }

    mpfr_clear (high_prec);

  /* 80 bit precision:
      - 3 limbs for 32-bit archs;
      - 2 limbs for 64-bit archs */
  mpfr_init2 (high_prec, 80);
  mpfr_set_d (high_prec, 1.0, MPFR_RNDD);

  hash_low_prec = mpfr_hash(low_prec);
  hash_high_prec = mpfr_hash(high_prec);

  if (hash_low_prec != hash_high_prec)
    {
      fprintf (stderr, "1.0 (prec. 1) and 1.0 (prec. 80) should have the "
                       "same hash.\n");
      exit (1);
    }

  mpfr_clear (low_prec);
  mpfr_clear (high_prec);

  double val = 1.0 / 3.0;

  mpfr_init2 (low_prec, 10);
  mpfr_set_d (low_prec, val, MPFR_RNDD);
  mpfr_init2 (high_prec, 50);
  mpfr_set_d (high_prec, val, MPFR_RNDD);

  hash_low_prec = mpfr_hash(low_prec);
  hash_high_prec = mpfr_hash(high_prec);

  if (hash_low_prec == hash_high_prec)
    {
      fprintf (stderr, "1.0 / 3.0 (prec.10) and 1.0 / 3.0 (prec. 50) should "
                       "not have the same hash");
      exit (1);
    }

  mpfr_clear (low_prec);
  mpfr_clear (high_prec);

  mpfr_free_cache ();
}

static void
test_constants (void)
{
  /* Arbitrary precision. */
  mpfr_prec_t p = 50;
  mpfr_t pi, log2, euler, catalan;
  mpfr_digest_t h_pi, h_log2, h_euler, h_catalan;

  mpfr_init2 (pi, p);
  mpfr_const_pi (pi, MPFR_RNDD);
  h_pi = mpfr_hash (pi);
  if (h_pi != H_PI)
    {
      fprintf (stderr, "pi digest should be " H_FORMAT "; got " H_FORMAT "\n",
               H_PI, h_pi);
      exit (1);
    }

  mpfr_init2 (log2, p);
  mpfr_const_log2 (log2, MPFR_RNDD);
  h_log2 = mpfr_hash (log2);
  if (h_log2 != H_LOG2)
    {
      fprintf (stderr, "log2 digest should be " H_FORMAT "; got " H_FORMAT "\n",
               H_LOG2, h_log2);
      exit (1);
    }

  mpfr_init2 (euler, p);
  mpfr_const_euler (euler, MPFR_RNDD);
  h_euler = mpfr_hash (euler);
  if (h_euler != H_EULER)
    {
      fprintf (stderr, "euler digest should be " H_FORMAT "; got " H_FORMAT "\n",
               H_EULER, h_euler);
      exit (1);
    }

  mpfr_init2 (catalan, p);
  mpfr_const_catalan (catalan, MPFR_RNDD);
  h_catalan = mpfr_hash (catalan);
  if (h_catalan != H_CATALAN)
    {
      fprintf (stderr, "catalan digest should be " H_FORMAT "; got " H_FORMAT "\n",
              H_CATALAN, h_catalan);
      exit (1);
    }

  mpfr_clear (pi);
  mpfr_clear (log2);
  mpfr_clear (euler);
  mpfr_clear (catalan);
  mpfr_free_cache ();
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  /* Special numbers */
  test_zero ();
  test_inf ();
  test_nan ();

  test_precision ();

  /* Those hard-coded digests of some mathematical constants. They guarantees
     that the hash function returns the same digest on both little endian and
     big endian machines, as long as they have the same exponent size. */
  test_constants ();

  tests_end_mpfr ();
}
