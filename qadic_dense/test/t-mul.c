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

    printf("mul... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing: a = a * b */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        qadic_dense_init(c);

        qadic_dense_randtest(a, state, ctx);
        qadic_dense_randtest(b, state, ctx);

        qadic_dense_mul(c, a, b, ctx);
        qadic_dense_mul(a, a, b, ctx);

        result = (qadic_dense_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_dense_print_pretty(c, ctx), printf("\n");
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);
        qadic_dense_clear(c);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    /* Check aliasing: b = a * b */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        qadic_dense_init(c);

        qadic_dense_randtest(a, state, ctx);
        qadic_dense_randtest(b, state, ctx);

        qadic_dense_mul(c, a, b, ctx);
        qadic_dense_mul(b, a, b, ctx);

        result = (qadic_dense_equal(b, c));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_dense_print_pretty(c, ctx), printf("\n");
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);
        qadic_dense_clear(c);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    /* Check aliasing: a = a + a */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(c);

        qadic_dense_randtest(a, state, ctx);

        qadic_dense_add(c, a, a, ctx);
        qadic_dense_add(a, a, a, ctx);

        result = (qadic_dense_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("c = "), qadic_dense_print_pretty(c, ctx), printf("\n");
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(c);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    /* Check that a * b == b * a */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, c1, c2;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        qadic_dense_init(c1);
        qadic_dense_init(c2);

        qadic_dense_randtest(a, state, ctx);
        qadic_dense_randtest(b, state, ctx);

        qadic_dense_mul(c1, a, b, ctx);
        qadic_dense_mul(c2, b, a, ctx);

        result = (qadic_dense_equal(c1, c2));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a  = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b  = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("c1 = "), qadic_dense_print_pretty(c1, ctx), printf("\n");
            printf("c2 = "), qadic_dense_print_pretty(c2, ctx), printf("\n");
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);
        qadic_dense_clear(c1);
        qadic_dense_clear(c2);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    /* Check that (a * b) * c == a * (b * c) for integral values */
    for (i = 0; i < 2000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, c, lhs, rhs;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        qadic_dense_init(c);
        qadic_dense_init(lhs);
        qadic_dense_init(rhs);

        qadic_dense_randtest_int(a, state, ctx);
        qadic_dense_randtest_int(b, state, ctx);
        qadic_dense_randtest_int(c, state, ctx);

        qadic_dense_mul(lhs, a, b, ctx);
        qadic_dense_mul(lhs, lhs, c, ctx);
        qadic_dense_mul(rhs, b, c, ctx);
        qadic_dense_mul(rhs, a, rhs, ctx);

        result = (qadic_dense_equal(lhs, rhs));
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a   = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b   = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("c   = "), qadic_dense_print_pretty(c, ctx), printf("\n");
            printf("lhs = "), qadic_dense_print_pretty(lhs, ctx), printf("\n");
            printf("rhs = "), qadic_dense_print_pretty(rhs, ctx), printf("\n");
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);
        qadic_dense_clear(c);
        qadic_dense_clear(lhs);
        qadic_dense_clear(rhs);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

