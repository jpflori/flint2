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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("set_ui_equal....");
    fflush(stdout);

    flint_randinit(state);

    /* equal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        ulong n;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        n = n_randtest(state);

        fmpz_poly_set_ui(a, n);
        fmpz_poly_set(b, a);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("n = %lu\n\n", n);
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        ulong m, n;

        fmpz_poly_init(a);
        fmpz_poly_init(b);

        m = n_randtest(state);
        n = n_randtest(state);
        while (m == n)
            n = n_randtest(state);
        fmpz_poly_set_ui(a, m);
        fmpz_poly_set_ui(b, n);

        result = (!fmpz_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("m = %lu\n\n", m);
            printf("n = %lu\n\n", n);
            printf("a = "), fmpz_poly_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
