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

    printf("mod_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f;
        ulong x, r1, r2;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(d, a);
        x = n_randtest_not_zero(state);

        r1 = fmpz_mod_ui(b, a, x);
        r2 = mpz_fdiv_r_ui(e, d, x);

        fmpz_get_mpz(f, b);

        result = ((mpz_cmp(e, f) == 0) && (r1 == r2));
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, x = %lu, r1 = %lu, r2 = %lu\n", d,
                 e, f, x, r1, r2);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
    }

    /* Check aliasing of a and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t d, e, f;
        ulong x, r1, r2;

        fmpz_init(a);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(d, a);
        x = n_randtest_not_zero(state);

        r1 = fmpz_mod_ui(a, a, x);
        r2 = mpz_fdiv_r_ui(e, d, x);

        fmpz_get_mpz(f, a);

        result = ((mpz_cmp(e, f) == 0) && (r1 == r2));
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf
                ("d = %Zd, e = %Zd, f = %Zd, x = %lu, r1 = %lu, r2 = %lu\n", d,
                 e, f, x, r1, r2);
            abort();
        }

        fmpz_clear(a);

        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
