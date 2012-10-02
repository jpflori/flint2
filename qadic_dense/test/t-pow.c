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

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "qadic_dense.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("pow... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing: a = a^e */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b;
        fmpz_t e;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        fmpz_init(e);

        qadic_dense_randtest(a, state, ctx);
        fmpz_randtest_unsigned(e, state, 6);

        qadic_dense_pow(b, a, e, ctx);
        qadic_dense_pow(a, a, e, ctx);

        result = (qadic_dense_equal(a, b));
        if (!result)
        {
            printf("FAIL (alias):\n\n");
            printf("a = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);
        fmpz_clear(e);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    /* Compare with multiplication, for integral values */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, c;
        fmpz_t e, f;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        qadic_dense_init(c);
        fmpz_init(f);
        fmpz_init(e);

        qadic_dense_randtest_int(a, state, ctx);
        fmpz_randtest_unsigned(e, state, 6);

        qadic_dense_pow(b, a, e, ctx);
        qadic_dense_one(c, ctx);
        for (fmpz_one(f); fmpz_cmp(f, e) <= 0; fmpz_add_ui(f, f, 1))
        {
            qadic_dense_mul(c, c, a, ctx);
        }

        result = (qadic_dense_equal(b, c));
        if (!result)
        {
            printf("FAIL (cmp with mul):\n\n");
            printf("a = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("e = "), fmpz_print(e), printf("\n");
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);
        qadic_dense_clear(c);
        fmpz_clear(e);
        fmpz_clear(f);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

