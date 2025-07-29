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
static const mpfr_digest_t H_PI      = 126454210;
static const mpfr_digest_t H_LOG2    = 1348199468;
static const mpfr_digest_t H_EULER   = 3879845459;
static const mpfr_digest_t H_CATALAN = 887850484;

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

  pZero_hash = mpfr_hash32 (pZero);
  nZero_hash = mpfr_hash32 (nZero);

  /* Zeros with different signs should be equal */
  if (pZero_hash != nZero_hash)
    {
      printf ("hash for +0.0 and -0.0 should be equal.\n"
              "H(+0.0) = %lu; H(-0.0) = %lu\n",
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
      if (mpfr_hash32 (pos[i]) != mpfr_hash32 (pos[j])
        && mpfr_hash32 (neg[i]) != mpfr_hash32 (neg[j]))
      {
        printf ("All zeros should be hashed the same, regardless their "
                "precision. i: %d j: %d\n", i, j);
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

  pInf_hash = mpfr_hash32 (pInf);
  nInf_hash = mpfr_hash32 (nInf);

  /* H(+Inf) and H(-Inf) should not be equal */
  if (pInf_hash == nInf_hash)
    {
      printf ("hash for +Inf and -Inf should not be equal.\n"
              "H(+Inf) = %lu; H(-Inf) = %lu\n",
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
    if (mpfr_hash32 (pos[i]) != mpfr_hash32 (pos[j]))
      {
        printf ("All +Inf should be hashed the same, regardless their "
                "precision. i: %d j: %d\n", i, j);
        exit (1);
      }
    if (mpfr_hash32 (neg[i]) != mpfr_hash32 (neg[j]))
      {
        printf ("All -Inf should be hashed the same, regardless their "
                "precision. i: %d j:%d\n", i, j);
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

  nan_hash = mpfr_hash32 (nan);
  unconventional_nan_hash = mpfr_hash32 (unconventional_nan);

  /* mpfr_hash32 should ignore the sign of NAN */
  if (nan_hash != unconventional_nan_hash)
    {
      printf ("hash for +NAN and -NAN should be equal.\n"
              "H(+NAN) = %lu; H(-NAN) = %lu\n",
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
      if (mpfr_hash32 (pos[i]) != mpfr_hash32 (pos[j]))
       {
          printf ("All NANs should be hashed the same, regardless their "
                  "precision. i: %d j: %d\n", i, j);
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

  hash_low_prec = mpfr_hash32 (low_prec);
  hash_high_prec = mpfr_hash32 (high_prec);

  if (hash_low_prec != hash_high_prec)
    {
      printf ("1.0 (prec. 1) and 1.0 (prec. 20) should have the same hash.\n");
      exit (1);
    }

  mpfr_clear (high_prec);

  /* 50 bit precision:
       - 2 limbs for 32-bit archs;
       - 1 limb for 64-bit archs */
  mpfr_init2 (high_prec, 50);
  mpfr_set_d (high_prec, 1.0, MPFR_RNDD);

  hash_low_prec = mpfr_hash32 (low_prec);
  hash_high_prec = mpfr_hash32 (high_prec);

  if (hash_low_prec != hash_high_prec)
    {
      printf ("1.0 (prec. 1) and 1.0 (prec. 50) should have the same hash.\n");
      exit (1);
    }

  mpfr_clear (high_prec);

  /* 80 bit precision:
      - 3 limbs for 32-bit archs;
      - 2 limbs for 64-bit archs */
  mpfr_init2 (high_prec, 80);
  mpfr_set_d (high_prec, 1.0, MPFR_RNDD);

  hash_low_prec = mpfr_hash32 (low_prec);
  hash_high_prec = mpfr_hash32 (high_prec);

  if (hash_low_prec != hash_high_prec)
    {
      printf ("1.0 (prec. 1) and 1.0 (prec. 80) should have the same hash.\n");
      exit (1);
    }

  mpfr_clear (low_prec);
  mpfr_clear (high_prec);

  double val = 1.0 / 3.0;

  mpfr_init2 (low_prec, 10);
  mpfr_set_d (low_prec, val, MPFR_RNDD);
  mpfr_init2 (high_prec, 50);
  mpfr_set_d (high_prec, val, MPFR_RNDD);

  hash_low_prec = mpfr_hash32 (low_prec);
  hash_high_prec = mpfr_hash32 (high_prec);

  if (hash_low_prec == hash_high_prec)
    {
      printf ("1.0 / 3.0 (prec. 10) and 1.0 / 3.0 (prec. 50) should "
              "not have the same hash.\n");
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
  h_pi = mpfr_hash32 (pi);
  if (h_pi != H_PI)
    {
      printf ("pi digest should be %lu; got %lu\n",
              H_PI, h_pi);
      exit (1);
    }

  mpfr_init2 (log2, p);
  mpfr_const_log2 (log2, MPFR_RNDD);
  h_log2 = mpfr_hash32 (log2);
  if (h_log2 != H_LOG2)
    {
      printf ("log2 digest should be %lu; got %lu\n",
              H_LOG2, h_log2);
      exit (1);
    }

  mpfr_init2 (euler, p);
  mpfr_const_euler (euler, MPFR_RNDD);
  h_euler = mpfr_hash32 (euler);
  if (h_euler != H_EULER)
    {
      printf ("euler digest should be %lu; got %lu\n",
              H_EULER, h_euler);
      exit (1);
    }

  mpfr_init2 (catalan, p);
  mpfr_const_catalan (catalan, MPFR_RNDD);
  h_catalan = mpfr_hash32 (catalan);
  if (h_catalan != H_CATALAN)
    {
      printf ("catalan digest should be %lu; got %lu\n",
              H_CATALAN, h_catalan);
      exit (1);
    }

  mpfr_clear (pi);
  mpfr_clear (log2);
  mpfr_clear (euler);
  mpfr_clear (catalan);
  mpfr_free_cache ();
}

static void
test_incremental_hashing (void)
{
  mpfr_digest_ctx_t ctx;
  const char *chunk1, *chunk2;
  size_t chunk1_len, chunk2_len;
  mpfr_digest_t got, expected;

  chunk1 = "Calculate the digest of ";
  chunk1_len = strlen (chunk1);
  chunk2 = "chunked bytes.";
  chunk2_len = strlen (chunk2);
  expected = 708666724;

  mpfr_digest_init (&ctx, MPFR_FNV_HASH_BYTES, mpfr_hash32_update,
                    mpfr_hash32_final);

  if (!mpfr_digest_update (&ctx, (const unsigned char *) chunk1, chunk1_len))
    {
      fprintf (stderr, "cannot calculate hash of chunk 1: \"%s\"\n", chunk1);
      exit (1);
    }

  if (!mpfr_digest_update (&ctx, (const unsigned char *) chunk2, chunk2_len))
    {
      fprintf (stderr, "cannot calculate hash of chunk 2: \"%s\"\n", chunk2);
      exit (1);
    }

  if (!mpfr_digest_final (&ctx, &got))
    {
      fprintf (stderr, "cannot get the resulting digest of chunk1 + chunk2:\n"
               "\"%s\" + \"%s\"\n", chunk1, chunk2);
      exit (1);
    }

  if (got != expected)
    {
      printf ("hash of chunk1 + chunk2 should be %lu, got %lu instead\n",
              expected, got);
      exit (1);
    }
}

static void
test_pi_incremental_hashing (void)
{
  mpfr_prec_t p = 50;
  mpfr_t pi;
  mpfr_digest_t h_pi;
  mpfr_digest_ctx_t ctx;

  mpfr_init2 (pi, p);
  mpfr_const_pi (pi, MPFR_RNDD);

  mpfr_digest_init (&ctx, MPFR_FNV_HASH_BYTES, mpfr_hash32_update,
                    mpfr_hash32_final);

  if (!mpfr_digest_update_m (&ctx, pi))
    {
      fprintf (stderr, "cannot calculate hash of pi constant with "
               "mpfr_digest_update_m\n");
      exit (1);
    }

  if (!mpfr_digest_final (&ctx, &h_pi))
    {
      fprintf (stderr, "cannot get the resulting digest of pi\n");
      exit (1);
    }

  if (h_pi != mpfr_hash32 (pi))
    {
      printf ("pi digest should be %lu; got %lu\n",
              H_PI, h_pi);
      exit (1);
    }

  mpfr_clear (pi);
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

  /* Incremental hashing of bytes sequences */
  test_incremental_hashing ();

  /* Incremental hashing of pi constant */
  test_pi_incremental_hashing ();

  tests_end_mpfr ();
}
