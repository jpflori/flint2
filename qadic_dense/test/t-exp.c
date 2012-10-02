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

#include "qadic_dense.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("exp... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, c;
        int ans1, ans2;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        qadic_dense_init(c);

        qadic_dense_randtest(a, state, ctx);
        qadic_dense_set(b, a);

        ans1 = qadic_dense_exp(c, b, ctx);
        ans2 = qadic_dense_exp(b, b, ctx);

        result = ((ans1 == ans2) && (!ans1 || qadic_dense_equal(b, c)));
        if (!result)
        {
            printf("FAIL (alias):\n\n");
            printf("a = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_dense_print_pretty(c, ctx), printf("\n");
            printf("ans1 = %d\n", ans1);
            printf("ans2 = %d\n", ans2);
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);
        qadic_dense_clear(c);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    /* Functional equation: exp(a + b) == exp(a) exp(b) */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long deg, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, c, d, e, f, g;
        int ans1, ans2, ans3;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        deg = n_randint(state, 10) + 1;
        N   = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, deg, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        qadic_dense_init(c);
        qadic_dense_init(d);
        qadic_dense_init(e);
        qadic_dense_init(f);
        qadic_dense_init(g);

        qadic_dense_randtest(a, state, ctx);
        qadic_dense_randtest(b, state, ctx);
        qadic_dense_add(c, a, b, ctx);

        ans1 = qadic_dense_exp(d, a, ctx);
        ans2 = qadic_dense_exp(e, b, ctx);
        qadic_dense_mul(f, d, e, ctx);

        ans3 = qadic_dense_exp(g, c, ctx);

        result = (!ans1 || !ans2 || (ans3 && qadic_dense_equal(f, g)));
        if (!result)
        {
            printf("FAIL (functional equation):\n\n");
            printf("a                 = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b                 = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("c = a + b         = "), qadic_dense_print_pretty(c, ctx), printf("\n");
            printf("d = exp(a)        = "), qadic_dense_print_pretty(d, ctx), printf("\n");
            printf("e = exp(b)        = "), qadic_dense_print_pretty(e, ctx), printf("\n");
            printf("f = exp(a) exp(b) = "), qadic_dense_print_pretty(f, ctx), printf("\n");
            printf("g = exp(a + b)    = "), qadic_dense_print_pretty(g, ctx), printf("\n");
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);
        qadic_dense_clear(c);
        qadic_dense_clear(d);
        qadic_dense_clear(e);
        qadic_dense_clear(f);
        qadic_dense_clear(g);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

