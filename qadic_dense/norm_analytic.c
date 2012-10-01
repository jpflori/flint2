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

/*
    Computes the norm of an element $x$ of $\mathbf{Z}_q$ via the identity 

        $\Norm(x) = \exp \Trace \log (x)$

    whenever $y = 1-x$ has valuation $v$ greater than $(p-1)^{-1}$.

    Assumes that $y$ is non-zero.
 */

void _qadic_dense_norm_analytic(fmpz_t rop, const fmpz *y, long v, long len, 
                          const fmpz *mod, const fmpz *invmod, long lenmod, 
                          const fmpz_t p, long N)
{
    const long d = lenmod - 1;
    fmpz_t pN, tru;
    long trv;
    fmpz *lg;

    fmpz_init(pN);
    fmpz_init(tru);
    lg = _fmpz_vec_init(d);

    fmpz_pow_ui(pN, p, N);

    _qadic_dense_log(lg, y, v, len, mod, invmod, lenmod, p, N, pN);

    _qadic_dense_trace(tru, lg, d, mod, invmod, lenmod, pN);

    if (!fmpz_is_zero(tru))
    {
        trv = fmpz_remove(tru, tru, p);
        _padic_exp(rop, tru, trv, p, N);
        fmpz_mod(rop, rop, pN);
    }
    else
    {
        fmpz_one(rop);
    }

    fmpz_clear(pN);
    fmpz_clear(tru);
    _fmpz_vec_clear(lg, d);
}

void qadic_dense_norm_analytic(padic_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx)
{
    const long N  = (&ctx->pctx)->N;
    const long d  = qadic_dense_ctx_degree(ctx);
    const fmpz *p = (&ctx->pctx)->p;

    /* N(p^v u) = p^{dv} N(u) */

    if (qadic_dense_is_zero(op) || d * op->val >= N)
    {
        padic_zero(rop);
    }
    else if (op->length == 1)
    {
        fmpz_t pN;
        int alloc;

        alloc = _padic_ctx_pow_ui(pN, N - d * op->val, (&ctx->pctx));

        fmpz_powm_ui(padic_unit(rop), op->coeffs + 0, d, pN);
        padic_val(rop) = d * op->val;

        if (alloc)
            fmpz_clear(pN);
    }
    else  /* len >= 2 */
    {
        fmpz *y;
        long w;

        y = _fmpz_vec_init(op->length);

        _fmpz_vec_neg(y, op->coeffs, op->length);
        fmpz_add_ui(y + 0, y + 0, 1);
        w = _fmpz_vec_ord_p(y, op->length, p);

        if ((w < 2 && *p == 2L) || w < 1)
        {
            printf("ERROR (qadic_dense_norm_analytic).  w = %ld.\n", w);
            abort();
        }

        _qadic_dense_norm_analytic(padic_unit(rop), y, w, op->length,
                                   ctx->mod->coeffs, ctx->invmod->coeffs, d + 1,
                                   p, N - d * op->val);
        padic_val(rop) = d * op->val;

        _fmpz_vec_clear(y, op->length);
    }
}

