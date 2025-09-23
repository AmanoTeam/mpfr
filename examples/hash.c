/* Hashing functions usage.

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

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>

int
main (int argc, char *argv[])
{
  mpfr_prec_t p = 50;
  mpfr_t pi;
  mpfr_digest_t h_pi;
  mpfr_digest_ctx_t ctx;

  mpfr_init2 (pi, p);
  mpfr_const_pi (pi, MPFR_RNDD);

  mpfr_digest_init (&ctx, MPFR_FNV_HASH32_BYTES, mpfr_hash32_update,
                    mpfr_hash32_final);

  if (!mpfr_digest_update_m (&ctx, pi))
    {
      fprintf (stderr, "An error occurred while executing digest update\n");
      goto cleanup;
    }

  if (!mpfr_digest_final (&ctx, &h_pi))
    {
      fprintf (stderr, "cannot get the resulting digest of pi\n");
      goto cleanup;
    }

  fprintf (stdout, "pi hash digest: %2lu\n", h_pi);

cleanup:
  mpfr_clear (pi);
  mpfr_free_cache ();
}
