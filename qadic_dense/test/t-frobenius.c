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

    printf("frobenius... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, c;
        long e;

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
        e = n_randint(state, 10) % d;

        qadic_dense_frobenius(c, b, e, ctx);
        qadic_dense_frobenius(b, b, e, ctx);

        result = (qadic_dense_equal(b, c));
        if (!result)
        {
            printf("FAIL (alias):\n\n");
            printf("a = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_dense_print_pretty(c, ctx), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);
        qadic_dense_clear(c);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    /* Check sigma^e(x) == x^{p^e} mod p for integral values */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, c, lhs, rhs;
        long e;

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
        e = n_randint(state, 10) % d;

        qadic_dense_frobenius(b, a, e, ctx);
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_pow_ui(t, p, e);
            qadic_dense_pow(c, a, t, ctx);
            fmpz_clear(t);
        }
        qadic_dense_scalar_mod_ppow(lhs, b, 1, ctx);
        qadic_dense_scalar_mod_ppow(rhs, c, 1, ctx);

        result = (qadic_dense_equal(lhs, rhs));
        if (!result)
        {
            printf("FAIL (sigma^e(x) = x^{p^e} mod p):\n\n");
            printf("a = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("c = "), qadic_dense_print_pretty(c, ctx), printf("\n");
            printf("lhs = "), qadic_dense_print_pretty(lhs, ctx), printf("\n");
            printf("rhs = "), qadic_dense_print_pretty(rhs, ctx), printf("\n");
            printf("e = %ld\n", e);
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

    /* Check sigma^e(x + y) = sigma^e(x) + sigma^e(y) on Zq */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, s, s1, s2, lhs, rhs;
        long e;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 10) + 1;
        qadic_dense_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        qadic_dense_init(s);
        qadic_dense_init(s1);
        qadic_dense_init(s2);
        qadic_dense_init(lhs);
        qadic_dense_init(rhs);

        qadic_dense_randtest_int(a, state, ctx);
        qadic_dense_randtest_int(b, state, ctx);
        e = n_randint(state, 10) % d;

        qadic_dense_add(s, a, b, ctx);
        qadic_dense_frobenius(lhs, s, e, ctx);
        qadic_dense_frobenius(s1, a, e, ctx);
        qadic_dense_frobenius(s2, b, e, ctx);
        qadic_dense_add(rhs, s1, s2, ctx);

        result = (qadic_dense_equal(lhs, rhs));
        if (!result)
        {
            printf("FAIL (sigma(a+b) = sigma(a) + sigma(b)):\n\n");
            printf("a = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("s = "), qadic_dense_print_pretty(s, ctx), printf("\n");
            printf("s1 = "), qadic_dense_print_pretty(s1, ctx), printf("\n");
            printf("s2 = "), qadic_dense_print_pretty(s2, ctx), printf("\n");
            printf("lhs = "), qadic_dense_print_pretty(lhs, ctx), printf("\n");
            printf("rhs = "), qadic_dense_print_pretty(rhs, ctx), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);
        qadic_dense_clear(s);
        qadic_dense_clear(s1);
        qadic_dense_clear(s2);
        qadic_dense_clear(lhs);
        qadic_dense_clear(rhs);

        fmpz_clear(p);
        qadic_dense_ctx_clear(ctx);
    }

    /* Check sigma^e(x * y) = sigma^e(x) * sigma^e(y) on Zq */
    for (i = 0; i < 1000; i++)
    {
        fmpz_t p;
        long d, N;
        qadic_dense_ctx_t ctx;

        qadic_dense_t a, b, s, s1, s2, lhs, rhs;
        long e;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 10) + 1;
        qadic_dense_ctx_init_conway(ctx, p, d, N, "a", PADIC_SERIES);

        qadic_dense_init(a);
        qadic_dense_init(b);
        qadic_dense_init(s);
        qadic_dense_init(s1);
        qadic_dense_init(s2);
        qadic_dense_init(lhs);
        qadic_dense_init(rhs);

        qadic_dense_randtest_int(a, state, ctx);
        qadic_dense_randtest_int(b, state, ctx);
        e = n_randint(state, 10) % d;

        qadic_dense_mul(s, a, b, ctx);
        qadic_dense_frobenius(lhs, s, e, ctx);
        qadic_dense_frobenius(s1, a, e, ctx);
        qadic_dense_frobenius(s2, b, e, ctx);
        qadic_dense_mul(rhs, s1, s2, ctx);

        result = (qadic_dense_equal(lhs, rhs));
        if (!result)
        {
            printf("FAIL (sigma(a*b) = sigma(a) * sigma(b)):\n\n");
            printf("a = "), qadic_dense_print_pretty(a, ctx), printf("\n");
            printf("b = "), qadic_dense_print_pretty(b, ctx), printf("\n");
            printf("s = "), qadic_dense_print_pretty(s, ctx), printf("\n");
            printf("s1 = "), qadic_dense_print_pretty(s1, ctx), printf("\n");
            printf("s2 = "), qadic_dense_print_pretty(s2, ctx), printf("\n");
            printf("lhs = "), qadic_dense_print_pretty(lhs, ctx), printf("\n");
            printf("rhs = "), qadic_dense_print_pretty(rhs, ctx), printf("\n");
            printf("e = %ld\n", e);
            abort();
        }

        qadic_dense_clear(a);
        qadic_dense_clear(b);
        qadic_dense_clear(s);
        qadic_dense_clear(s1);
        qadic_dense_clear(s2);
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

