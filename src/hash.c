/* mpfr_hash32 -- function to calculate 32-bit hash digest for MPFR numbers

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

#include "mpfr-impl.h"

#define FNV32_PRIME 0x01000193U

/* To extract bytes from a limb we need to check for endianness */
#if defined (HAVE_LITTLE_ENDIAN)
#define limb_byte_index(b) (sizeof (mp_limb_t) - 1 - (b))
#elif defined (HAVE_BIG_ENDIAN)
#define limb_byte_index(b) (b)
#endif

/* Those constants contain the three bytes for the 4 special numbers +/-0,
   NAN, +Inf, -Inf, stored in little endian. The first byte represents the
   sign (0 for positive and 1 for negative), the second and the third bytes
   contain the 16-bit little endian exponent. Neither precision nor
   significand are encoded for these numbers, because mpfr_min_prec (x) = 0
   where x in {+/-0,NAN,+Inf,-Inf} */
#define MPFR_SINGULAR_DIGEST_SIZE 3
static const unsigned char default_zero[]    = { 0x00, 0x01, 0x80 };
static const unsigned char default_nan[]     = { 0x00, 0x02, 0x80 };
static const unsigned char default_pos_inf[] = { 0x00, 0x03, 0x80 };
static const unsigned char default_neg_inf[] = { 0x01, 0x03, 0x80 };

/* Removing the leading zeros of a number guarantees that the encoding will
   be the same on different architectures with different sizes of the
   numbers */
static size_t
count_relevant_bytes (void *x, size_t size)
{
  unsigned char *bytes = (unsigned char *)x;
#if defined (HAVE_BIG_ENDIAN)
  for (; size > 1 && *bytes == 0; bytes++, size--) ;
  return size;
#elif defined (HAVE_LITTLE_ENDIAN)
  for (; size > 0 && bytes[--size] == 0;) ;
  return ++size;
#endif
}

/* copy in {bytes,n} the relevant bytes from x, and return the length n */
static size_t
le_relevant_bytes (unsigned char *bytes, size_t size, void *x)
{
  size_t relevant_size = count_relevant_bytes (x, size);
#if defined (HAVE_BIG_ENDIAN)
  unsigned char *src = (unsigned char *)x;
  src += (size - relevant_size);
  for (size_t i = 0; i < relevant_size; ++i)
    bytes[i] = src[relevant_size - 1 - i];
#elif defined (HAVE_LITTLE_ENDIAN)
  memcpy (bytes, x, relevant_size);
#endif
  return relevant_size;
}

/* return the number of bytes occupied by 'bits' bits, i.e., ceil(bits/8) */
static size_t
bytes_for (mpfr_prec_t bits)
{
  return (bits + 7) / 8;
}

/* return the number of bytes used by x */
static size_t
get_bytes_size (mpfr_srcptr x)
{
  return 1                                            /* sign size */
         + sizeof (mpfr_prec_t)                       /* precision size */
         + sizeof (mpfr_exp_t)                        /* exponent size */
         + (MPFR_LIMB_SIZE (x) * sizeof (mp_limb_t)); /* significand size */
}

static digest32_t
fnv32 (digest32_t hash, const unsigned char *bytes, size_t bytes_len)
{
  size_t i;

  /* if bytes = NULL, we should have bytes_len = 0 */

  for (i = 0; i < bytes_len; i++)
    {
      hash ^= bytes[i];
      hash *= FNV32_PRIME;
    }

  return hash;
}

static const unsigned char *
get_singular_number (mpfr_srcptr x)
{
  if (MPFR_IS_ZERO (x))
    return default_zero;

  if (MPFR_IS_NAN (x))
    return default_nan;

  if (MPFR_IS_POS (x))
    return default_pos_inf;

  MPFR_ASSERTD (MPFR_IS_NEG (x));
  return default_neg_inf;
}

/* put in bytes the unique bytes of a non-singular number */
static int
non_singular_unique_bytes (mpfr_srcptr x, mpfr_bytes_t *bytes)
{
  int i, j;
  size_t written_bytes, limb_bytes = 0, min_prec_byte_size = 0,
         limb_size = 0, bytes_size = 0;
  unsigned char *mpfr_bytes = NULL, *l_bytes = NULL;
  unsigned char sign;
  mpfr_exp_t exp;
  mpfr_prec_t prec;
  mp_limb_t *limbs = NULL;

  sign = MPFR_IS_NEG (x) ? 1 : 0;
  exp = MPFR_GET_EXP (x);

  /* We only encode the minimum number of bits required to
     represent the number */
  prec = mpfr_min_prec (x);
  /* since x is non-singular, we should have prec >= 1 */
  MPFR_ASSERTD (prec > 0);
  /* The minimum number of bytes required to represent prec bits */
  min_prec_byte_size = bytes_for (prec);
  limb_size = MPFR_LIMB_SIZE (x);
  limbs = MPFR_MANT (x);

  /* The size of the byte array required to encode the number is
     always less (or equal) than `bytes_size`:
       min_prec_byte_size <= bytes_size */
  bytes_size = get_bytes_size (x);
  /* again, since x is non-singular, we have bytes_size > 0 */
  MPFR_ASSERTD (bytes_size > 0);
  mpfr_bytes = (unsigned char *) malloc (bytes_size);

  /* We encode (sequentially):
       - the sign (written_bytes = 1);
       - the precision (written_bytes <= sizeof (mpfr_prec_t));
       - the exponent (written_bytes <= sizeof (mpfr_exp_t)) */
  written_bytes = le_relevant_bytes (mpfr_bytes, 1, &sign);
  written_bytes += le_relevant_bytes (mpfr_bytes + written_bytes,
                                      sizeof (mpfr_prec_t), &prec);
  written_bytes += le_relevant_bytes (mpfr_bytes + written_bytes,
                                      sizeof (mpfr_exp_t), &exp);

  /* To encode the significand, we need to encode all the limbs in a reverse
     order. If the precision is not a multiple of the size of `mp_limb_t`,
     this algorithm should get rid of all the leading zeros. It should copy
     limbs' bytes from the least to the most significant. */
  for (i = limb_size - 1; i >= 0 && limb_bytes < min_prec_byte_size; --i)
    {
      l_bytes = (unsigned char *) (limbs + i);

      for (j = 0; j < sizeof (mp_limb_t) && limb_bytes < min_prec_byte_size; ++j)
        {
          mpfr_bytes[written_bytes++] = l_bytes[limb_byte_index (j)];
          limb_bytes++;
        }
    }
  bytes->content = mpfr_bytes;
  bytes->len = written_bytes;

  return 1;
}

/* put in bytes the unique bytes of x */
int
mpfr_unique_bytes (mpfr_srcptr x, mpfr_bytes_t *bytes)
{
  const unsigned char *singular_num;

  if (x == NULL || bytes == NULL)
    return 0;

  if (MPFR_IS_SINGULAR (x))
    {
      bytes->content = (unsigned char *) malloc (MPFR_SINGULAR_DIGEST_SIZE);
      if (!bytes->content)
      {
        bytes->len = 0;
        return 0;
      }
      bytes->len = MPFR_SINGULAR_DIGEST_SIZE;

      singular_num = get_singular_number (x);
      memcpy (bytes->content, singular_num, MPFR_SINGULAR_DIGEST_SIZE);
      return 1;
    }

  return non_singular_unique_bytes (x, bytes);
}

/* free the bytes array */
void
mpfr_bytes_free (mpfr_bytes_t *bytes)
{
  free (bytes->content);
  bytes->content = NULL;
  bytes->len = 0;
}

/* main hash function */
int
mpfr_hash32 (mpfr_digest_t *digest, mpfr_srcptr x)
{
  mpfr_bytes_t bytes = { 0 };

  if (!mpfr_unique_bytes (x, &bytes))
    return 0;

  *digest = fnv32 (MPFR_HASH32_BASIS, bytes.content, bytes.len);
  mpfr_bytes_free (&bytes);

  return 1;
}

int
mpfr_hash32_update (mpfr_digest_ctx_t ctx, const unsigned char *bytes,
                    size_t l)
{
  if (!ctx || !bytes)
    return 0;

  ctx->hash = fnv32 (ctx->hash, bytes, l);

  return 1;
}

int
mpfr_hash32_final (const mpfr_digest_ctx_t ctx, mpfr_digest_t *digest)
{
  *digest = ctx->hash;
  return 1;
}

mpfr_digest_ctx_t
mpfr_digest_init (mpfr_digest_t basis, size_t digest_size,
                  hash_update_fn_t update_fn, hash_final_fn_t final_fn)
{
  mpfr_digest_ctx_t ctx;

  ctx = (mpfr_digest_ctx_t) mpfr_allocate_func (sizeof (*ctx));
  ctx->hash = basis;
  ctx->digest_size = digest_size;
  ctx->update_fn = update_fn;
  ctx->final_fn = final_fn;

  return ctx;
}

void
mpfr_digest_ctx_clear (mpfr_digest_ctx_t *ctx)
{
  mpfr_free_func (*ctx, sizeof (**ctx));
  *ctx = NULL;
}

int
mpfr_digest_update (mpfr_digest_ctx_t ctx, const unsigned char *data,
                    size_t len)
{
  return ctx->update_fn (ctx, data, len);
}

int
mpfr_digest_update_m (mpfr_digest_ctx_t ctx, mpfr_srcptr x)
{
  int ret = 1;
  mpfr_bytes_t bytes;

  if (!mpfr_unique_bytes (x, &bytes))
    return 0;

  ret = ctx->update_fn (ctx, bytes.content, bytes.len);
  mpfr_bytes_free (&bytes);

  return ret;
}

int
mpfr_digest_final (const mpfr_digest_ctx_t ctx, mpfr_digest_t *digest)
{
  return ctx->final_fn (ctx, digest);
}

mpfr_digest_t
mpfr_digest_ctx_get_hash (const mpfr_digest_ctx_t ctx)
{
  return ctx->hash;
}

void
mpfr_digest_ctx_set_hash (mpfr_digest_ctx_t ctx, mpfr_digest_t hash)
{
  ctx->hash = hash;
}
