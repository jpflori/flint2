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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("mulhigh_classical....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        long len, start;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);

        len = b->length + c->length - 1;
        start = (len <= 0) ? 0 : n_randint(state, b->length + c->length);

        fmpz_poly_mulhigh_classical(a, b, c, start);
        fmpz_poly_mulhigh_classical(b, b, c, start);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        long len, start;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);

        len = b->length + c->length - 1;
        start = (len <= 0) ? 0 : n_randint(state, b->length + c->length - 1);

        fmpz_poly_mulhigh_classical(a, b, c, start);
        fmpz_poly_mulhigh_classical(c, b, c, start);

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(c), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Compare with mul_basecase */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d;
        long len, start;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);

        len = b->length + c->length - 1;
        start = (len <= 0) ? 0 : n_randint(state, b->length + c->length - 1);

        fmpz_poly_mul_classical(a, b, c);
        if (a->length >= start)
            _fmpz_vec_zero(a->coeffs, start);
        fmpz_poly_mulhigh_classical(d, b, c, start);

        result = (fmpz_poly_equal(a, d));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_print(a), printf("\n\n");
            fmpz_poly_print(d), printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
