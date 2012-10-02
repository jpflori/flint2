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

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>

#include "qadic_dense.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("log_rectangular... ");
    fflush(stdout);

    flint_randinit(state);

/** p == 2 *******************************************************************/

/** p > 2 ********************************************************************/

    /* Check aliasing: a = log(a) */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        long ans1, ans2;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);

        qadic_dense_randtest_not_zero(a, state, ctx);
        if (a->val < 1)
            a->val = 1;
        padic_poly_reduce(a, &ctx->pctx);
        qadic_dense_one(b, ctx);
        qadic_dense_add(a, a, b, ctx);

        ans1 = qadic_dense_log_rectangular(b, a, ctx);
        ans2 = qadic_dense_log_rectangular(a, a, ctx);

        result = (ans1 == ans2) && (!ans1 || qadic_dense_equal(a, b));
        if (!result)
        {
            printf("FAIL (aliasing):\n\n");
            printf("a = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    /* Check: log(a) + log(b) == log(a * b) */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long deg, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, c, d, e, f, g;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 3 + n_randint(state, 3), 1));
        deg = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, deg, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        qadic_dense_init(c);
        qadic_dense_init(d);
        qadic_dense_init(e);
        qadic_dense_init(f);
        qadic_dense_init(g);

        qadic_dense_randtest_not_zero(a, state, ctx);
        if (a->val < 1) 
            a->val = 1;
        padic_poly_reduce(a, &ctx->pctx);
        qadic_dense_one(c, ctx);
        qadic_dense_add(a, a, c, ctx);

        qadic_dense_randtest_not_zero(b, state, ctx);
        if (b->val < 1) 
            b->val = 1;
        padic_poly_reduce(b, &ctx->pctx);
        qadic_dense_one(c, ctx);
        qadic_dense_add(b, b, c, ctx);

        qadic_dense_mul(c, a, b, ctx);

        qadic_dense_log_rectangular(d, a, ctx);
        qadic_dense_log_rectangular(e, b, ctx);
        qadic_dense_add(f, d, e, ctx);

        qadic_dense_log_rectangular(g, c, ctx);

        result = (qadic_dense_equal(f, g));
        if (!result)
        {
            printf("FAIL (functional equation):\n\n");
            printf("a                   = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b                   = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("c = a * b           = "), qadic_dense_print_pretty(c, ctx), printf("\n");
            printf("d = log(a)          = "), qadic_dense_print_pretty(d, ctx), printf("\n");
            printf("e = log(b)          = "), qadic_dense_print_pretty(e, ctx), printf("\n");
            printf("f = log(a) + log(b) = "), qadic_dense_print_pretty(f, ctx), printf("\n");
            printf("g = log(a * b)      = "), qadic_dense_print_pretty(g, ctx), printf("\n");
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

    /* Check: log(exp(x)) == x */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long deg, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, c;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        deg = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_dense_ctx_init_conway(ctx, p, deg, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        qadic_dense_init(c);

        qadic_dense_randtest(a, state, ctx);
        if (!qadic_dense_is_zero(a) && (a->val < 1 || (*p == 2L && a->val < 2)))
        {
            a->val = (*p == 2L) + 1;
            qadic_dense_scalar_mod_ppow(a, a, N, ctx);
        }

        qadic_dense_exp(b, a, ctx);
        qadic_dense_log_rectangular(c, b, ctx);

        result = (qadic_dense_equal(a, c));
        if (!result)
        {
            printf("FAIL (log(exp(x)) == x):\n\n");
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

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

