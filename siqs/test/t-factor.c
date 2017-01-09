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

    Copyright (C) 2015 Nitin Kumar

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "siqs.h"

void randprime(fmpz_t p, flint_rand_t state, slong bits)
{
    fmpz_randbits(p, state, bits);
 
    if (fmpz_sgn(p) < 0)
       fmpz_neg(p, p);

    if (fmpz_is_even(p))
       fmpz_add_ui(p, p, 1);
 
    while (!fmpz_is_probabprime(p))
       fmpz_add_ui(p, p, 2);
}

int main(void)
{
   slong i;
   fmpz_t n, x, y, z;
   fmpz_factor_t factors;
   FLINT_TEST_INIT(state);

   fmpz_init(x);
   fmpz_init(y);
   fmpz_init(z);
   fmpz_init(n);

   flint_printf("factor....");
   fflush(stdout);

   for (i = 0; i < 20; i++) /* Test random n, two factors */
   {
      randprime(x, state, 50);
      randprime(y, state, 50);

      fmpz_mul(n, x, y);

      fmpz_factor_init(factors);

      qsieve_factor(n, factors);

      if (factors->num < 2)
      {
         flint_printf("FAIL:\n");
         flint_printf("%ld factors found\n", factors->num);
         abort();
      }

      fmpz_factor_clear(factors);
   }

   for (i = 0; i < 20; i++) /* Test random n, three factors */
   {
      randprime(x, state, 40);
      randprime(y, state, 40);
      randprime(z, state, 40);

      fmpz_mul(n, x, y);
      fmpz_mul(n, n, z);

      fmpz_factor_init(factors);

      qsieve_factor(n, factors);

      if (factors->num < 3)
      {
         flint_printf("FAIL:\n");
         flint_printf("%ld factors found\n", factors->num);
         abort();
      }

      fmpz_factor_clear(factors);
   }

   for (i = 0; i < 20; i++) /* Test random n, small factors */
   {
      randprime(x, state, 10);
      randprime(y, state, 10);
      randprime(z, state, 40);

      fmpz_mul(n, x, y);
      fmpz_mul(n, n, z);

      fmpz_factor_init(factors);

      qsieve_factor(n, factors);

      if (factors->num < 3)
      {
         flint_printf("FAIL:\n");
         flint_printf("%ld factors found\n", factors->num);
         abort();
      }

      fmpz_factor_clear(factors);
   }

   for (i = 0; i < 5; i++) /* Test random squares */
   {
      randprime(x, state, 40);
      
      fmpz_mul(n, x, x);

      fmpz_factor_init(factors);

      qsieve_factor(n, factors);

      if (factors->num < 1)
      {
         flint_printf("FAIL:\n");
         flint_printf("%ld factors found\n", factors->num);
         abort();
      }

      fmpz_factor_clear(factors);
   }

   fmpz_clear(n);
   fmpz_clear(x);
   fmpz_clear(y);
   fmpz_clear(z);

   FLINT_TEST_CLEANUP(state);

   flint_printf("PASS\n");
   return 0;
}