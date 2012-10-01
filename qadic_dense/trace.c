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

#include "qadic_dense.h"

void _qadic_dense_trace(fmpz_t rop, const fmpz *op, long len, 
                  const fmpz *mod, const fmpz *invmod, long lenmod, const fmpz_t pN)
{
    const long d = lenmod - 1;

    long i, l;
    fmpz *t;

    t = _fmpz_vec_init(d);

    fmpz_set_ui(t + 0, d);
    for (i = 1; i < d; i++)
    {
        for (l = lenmod - 2; l >= 0 && l >= d - (i - 1); l--)
        {
            fmpz_addmul(t + i, t + (l + i - d), mod + l);
        }

        if (l >= 0 && l == d - i)
        {
            fmpz_addmul_ui(t + i, mod + l, i);
        }

        fmpz_neg(t + i, t + i);
        fmpz_mod(t + i, t + i, pN);
    }

    fmpz_zero(rop);
    for (i = 0; i < d; i++)
    {
        fmpz_addmul(rop, op + i, t + i);
    }
    fmpz_mod(rop, rop, pN);

    _fmpz_vec_clear(t, d);
}

void qadic_dense_trace(padic_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (qadic_dense_is_zero(op) || op->val >= N)
    {
        padic_zero(rop);
    }
    else
    {
        const int d = qadic_dense_ctx_degree(ctx);
        fmpz_t pN;
        int alloc;

        alloc = _padic_ctx_pow_ui(pN, N - op->val, &ctx->pctx);

        _qadic_dense_trace(padic_unit(rop), op->coeffs, op->length, 
                           ctx->mod->coeffs, ctx->invmod->coeffs, d + 1, pN);
        padic_val(rop) = op->val;

        _padic_canonicalise(rop, &ctx->pctx);

        if (alloc)
            fmpz_clear(pN);
    }
}

