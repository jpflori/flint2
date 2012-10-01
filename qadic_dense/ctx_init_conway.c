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

    Copyright (C) 2011 Sebastian Pancratz
 
******************************************************************************/

#include <stdio.h>
#include <string.h>

#include "fmpz_vec.h"
#include "padic.h"
#include "qadic_dense.h"

#ifndef FLINT_CPIMPORT
#define FLINT_CPIMPORT "/home/user/FLINT/flint-2/qadic_dense/CPimport.txt"
#endif

int fmpz_poly_init_conway(fmpz_poly_t poly,
                          const fmpz_t p, long d)
{
    char *buf;
    FILE *file;

    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        printf("Exception (fmpz_poly_init_conway).  Conway polynomials \n");
        printf("are only available for primes up to 109987.\n");
        return 1;
    }

    buf  = flint_malloc(832);
    file = fopen(FLINT_CPIMPORT, "r");

    if (!file)
    {
        printf("Exception (fmpz_poly_init_conway).  File loading.\n");
        return 2;
    }

    while (fgets(buf, 832, file))
    {
        char *tmp = buf;

        /* Different prime? */
        if (fmpz_cmp_ui(p, atoi(tmp)))
            continue;

        while (*tmp++ != ' ') ;

        /* Same degree? */
        if (d == atoi(tmp))
        {
            long i;
            char *ptr;

            /* Initialization */
            fmpz_poly_init2(poly, d + 1);

            /* Read coefficients */
            ptr = tmp;

            for (i = 0; i <= d; i++)
            {
                while (*ptr++ != ' ') ;

                fmpz_poly_set_coeff_si(poly, i, atoi(ptr));
            }

            fclose(file);
            flint_free(buf);
            return 0;
        }
    }

    fclose(file);
    flint_free(buf);

    printf("Exception (fmpz_poly_init_conway).  The polynomial for \n");
    printf("(p,d) = (%ld,%ld) is not present in the database.\n", *p, d);
    return 1;
}

int padic_poly_init_conway(padic_poly_t poly,
                           const fmpz_t p, long d,
                           const padic_ctx_t ctx)
{
    fmpz_poly_t cpoly;
    if (fmpz_poly_init_conway(cpoly, p, d))
    {
        return 1;
    }
    padic_poly_init2(poly, d + 1);
    padic_poly_set_fmpz_poly(poly, cpoly, ctx);
    fmpz_poly_clear(cpoly);
    return 0;
}

void qadic_dense_ctx_init_conway(qadic_dense_ctx_t ctx,
                           const fmpz_t p, long d, long N, const char *var, 
                           enum padic_print_mode mode)
{
    padic_ctx_init(&ctx->pctx, p, N, mode);

    if (padic_poly_init_conway(ctx->mod, p, d, &ctx->pctx)) {
        printf("Error.  Conway polynomial not found.\n");
        abort();
    }

    /* Precomputed inverse for fast modular reduction */
    qadic_dense_ctx_init_inv(ctx, &ctx->pctx, d, N);

    ctx->var = flint_malloc(strlen(var));
    strcpy(ctx->var, var);
    return;
}
