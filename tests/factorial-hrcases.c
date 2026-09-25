/* Find long sequences of 0s or 1s in n! (for tfactorial.c)

Copyright 2026 Free Software Foundation, Inc.
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

/* Output format: a list of C arrays whose elements are
 *   [0] The value of n.
 *   [1] The repeated bit in the sequence of 0s or 1s.
 *   [2] The length of this sequence.
 *   [3] The position of the MSB of this sequence (0 for the MSB of n!).
 *   [4] The size of n!, i.e. its exponent (useful for GNU MPFR).
 *
 * Note: stdout is flushed as in practice, this program will not output
 * many data. So this is very useful when stdout is piped, e.g. to tee.
 * This avoids the need to run this program via stdbuf.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>

int main (int argc, char **argv)
{
  unsigned long minlen, pmax;
  unsigned long n0, n1, n;
  char *end;
  mpz_t f;

  if (argc  < 5)
    {
      fprintf (stderr, "Usage: %s <minlen> <pmax> <n0> <n1>\n", argv[0]);
      exit (1);
    }

  minlen = strtoul (argv[1], &end, 10);
  if (*end != '\0')
    {
      fprintf (stderr, "%s: error on minlen\n", argv[0]);
      exit (1);
    }

  pmax = strtoul (argv[2], &end, 10);
  if (*end != '\0')
    {
      fprintf (stderr, "%s: error on pmax\n", argv[0]);
      exit (1);
    }

  n0 = strtoul (argv[3], &end, 10);
  if (*end != '\0')
    {
      fprintf (stderr, "%s: error on n0\n", argv[0]);
      exit (1);
    }

  n1 = strtoul (argv[4], &end, 10);
  if (*end != '\0')
    {
      fprintf (stderr, "%s: error on n1\n", argv[0]);
      exit (1);
    }

  mpz_init (f);
  mpz_fac_ui (f, n0 - 1);

  for (n = n0; n <= n1; n++)
    {
      unsigned long s, i, j;
      int b;

      mpz_mul_ui (f, f, n); /* f = n! */
      s = mpz_sizeinbase (f, 2);

#ifdef VERBOSE
      printf ("\n");
      mpz_out_str (stdout, 2, f);
      printf (" (%lu)\n", s);
      fflush (stdout);
#endif

      /* Goal of pmax: to restrict the precision before the sequence. */
      i = s > pmax ? s - pmax - 1 : 0;
      b = mpz_tstbit (f, i);
      while (i > 0 && mpz_tstbit (f, i-1) == b)
        i--;

      for (; i < s; i = j, b = !b)
        {
          assert (mpz_tstbit (f, i) == b);
          j = (b ? mpz_scan0 : mpz_scan1) (f, i);
          assert (j > i);
          if (i == 0)
            continue;
          unsigned long len = j - i;
          if (len < minlen)
            continue;
          printf ("{ %lu, %d, %lu, %lu, %lu },\n", n, b, len, s - j, s);
          fflush (stdout);
        }
    }

  mpz_clear (f);
  return 0;
}
