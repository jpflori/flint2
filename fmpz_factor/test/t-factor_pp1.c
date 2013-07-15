/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "ulong_extras.h"

int main(void)
{
   int i, j, result;
   ulong count = 0UL;
   gmp_randstate_t st;
   flint_rand_t state;
   gmp_randinit_default(st);
   flint_randinit(state);

   printf("factor_pp1....");
   fflush(stdout);

   for (i = 0; i < 50 * flint_test_multiplier(); i++) /* Test random numbers */
   {
      mp_bitcnt_t bits;
      mpz_t m, n;
      fmpz_t n1, n2, r;

      mpz_init(n);
      mpz_init(m);
      fmpz_init(n1);
      fmpz_init(n2);
      fmpz_init(r);

      do {
         mpz_urandomb(n, st, n_randint(state, 128) + 2); 
      } while (mpz_cmp_ui(n, 2) < 0);
      do {
         mpz_urandomb(m, st, n_randint(state, 50) + 2); 
      } while (!mpz_probab_prime_p(m, 20));
      mpz_mul(n, n, m);

      fmpz_set_mpz(n1, n);
      bits = FLINT_MIN(fmpz_bits(n1), FLINT_BITS);

      for (j = 0; j < 20; j++)
      {
         fmpz_factor_pp1(n2, n1, 10000, 10000, n_randbits(state, bits - 2) + 3);
         if (fmpz_cmp_ui(n2, 1) > 0) break;
      }
      
      if (fmpz_cmp_ui(n2, 1) > 0)
      {
         count++;
         fmpz_mod(r, n1, n2);
         result = (fmpz_is_zero(r));
         if (!result)
         {
            printf("FAIL:\n");
            printf("n1 = ");
            fmpz_print(n1);
            printf(", n2 = ");
            fmpz_print(n2);
            printf("\n"); 
            fmpz_print(r); printf("\n");
            abort();
         }
      }

      fmpz_clear(n1);
      fmpz_clear(n2);
      fmpz_clear(r);
      mpz_clear(m);
      mpz_clear(n);
   }
   
   if (count < 49 * flint_test_multiplier())
   {
      printf("FAIL:\n");
      printf("Only %lu numbers factored\n", count); 
      abort();
   }

   flint_randclear(state);
   gmp_randclear(st);

   printf("PASS\n");
   return 0;
}
