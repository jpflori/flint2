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

    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("combit....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong j;
        fmpz_t a;
        mpz_t b, c;

        fmpz_init(a);
        mpz_init(b);
        mpz_init(c);

        fmpz_randtest(a, state, 2 * FLINT_BITS);
        fmpz_get_mpz(b, a);
        j = n_randint(state, 3 * FLINT_BITS);

        fmpz_combit(a, j);
        mpz_combit(b, j);
        fmpz_get_mpz(c, a);

        result = (mpz_cmp(b, c) == 0);

        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_print(a), printf("\n");
            gmp_printf("b = %Zd\n", b);
            printf("j = %ld\n", j);
            abort();
        }

        fmpz_clear(a);
        mpz_clear(c);
        mpz_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
