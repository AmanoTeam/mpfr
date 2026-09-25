/* MPFR hash interface usage.

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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <mpfr.h>

#define DJB2_BASIS 0x00001505

static uint32_t
djb2(uint32_t hash, const unsigned char *bytes, size_t bytes_len)
{
  for (size_t i = 0; i < bytes_len; i++)
      hash = ((hash << 5) + hash) + bytes[i];

  return hash;
}

static int
djb2_update (mpfr_digest_ctx_t ctx, const unsigned char *bytes,
             size_t l)
{
  mpfr_digest_t updated_hash;

  if (!ctx || !bytes)
    return 0;

  updated_hash = djb2 (mpfr_digest_ctx_get_hash (ctx), bytes, l);
  mpfr_digest_ctx_set_hash (ctx, updated_hash);

  return 1;
}

static int
djb2_final (const mpfr_digest_ctx_t ctx, mpfr_digest_t *digest)
{
  *digest = mpfr_digest_ctx_get_hash (ctx);
  return 1;
}

static void
custom_mpfr_hash32 (mpfr_t x)
{
  mpfr_digest_t hash;
  mpfr_digest_ctx_t ctx;

  ctx = mpfr_digest_init (DJB2_BASIS, MPFR_HASH32_BYTES,
                          djb2_update, djb2_final);

  if (!mpfr_digest_update_m (ctx, x))
    {
      fprintf (stderr, "[custom_mpfr_hash32] An error occurred while "
                       "executing digest update\n");
      goto cleanup;
    }

  if (!mpfr_digest_final (ctx, &hash))
    {
      fprintf (stderr, "[custom_mpfr_hash32] cannot get the digest for x\n");
      goto cleanup;
    }

  fprintf (stdout, "custom 32-bit hash digest:\t%2lu\n", hash);

cleanup:
  mpfr_digest_ctx_clear (&ctx);
}

static void
default_mpfr_hash32 (mpfr_t x)
{
  mpfr_digest_t hash;
  mpfr_digest_ctx_t ctx;

  ctx = mpfr_digest_init (MPFR_HASH32_BASIS, MPFR_HASH32_BYTES,
                          mpfr_hash32_update, mpfr_hash32_final);

  if (!mpfr_digest_update_m (ctx, x))
    {
      fprintf (stderr, "[default_mpfr_hash32] An error occurred while "
                       "executing digest update\n");
      goto cleanup;
    }

  if (!mpfr_digest_final (ctx, &hash))
    {
      fprintf (stderr, "[default_mpfr_hash32] cannot get the digest for x\n");
      goto cleanup;
    }

  fprintf (stdout, "default 32-bit hash digest:\t%2lu\n", hash);

cleanup:
  mpfr_digest_ctx_clear (&ctx);
}

static void
hash32_shortcut (mpfr_t x)
{
  int ret;
  mpfr_digest_t hash;

  ret = mpfr_hash32 (&hash, x);

  if (!ret)
    {
      fprintf (stderr, "[hash32_shortcut] An error occurred while "
                       "executing digest update\n");
      return;
    }

  fprintf (stdout, "mpfr_hash32 digest (shortcut):\t%2lu\n", hash);
}

int
main (void)
{
  mpfr_prec_t p = 50;
  mpfr_t pi;

  mpfr_init2 (pi, p);
  mpfr_const_pi (pi, MPFR_RNDD);

  /* this function calculates the 32 bit hash of the number pi using the
     default algorithm implemented in MPFR. If you need to calculate a ready
     to use digest of a single number, this is the simplest way to do that */
  hash32_shortcut (pi);

  /* this is the standard way to compute hashed in MPFR. See mpfr.texi for
     further details on the APIs usage */
  default_mpfr_hash32 (pi);

  /* the same APIs can be used to compute the digest by using user-implemented
     algorithms. It uses a very simple implementation of the djb2 hash function
     (see the djb2 function here). As documented in mpfr.texi, you can implement
     your own hash function and pass it to MPFR as a callback to
     mpfr_digest_init */
  custom_mpfr_hash32 (pi);

cleanup:
  mpfr_clear (pi);
  mpfr_free_cache ();
}
