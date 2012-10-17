/*============================================================================

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
    Forms the product of (op1,len1) and (op2,len2) modulo (a,j,lena) and pN.
    Assumes that len1 >= len2 > 0.  Requires rop to be of size at least 
    len1 + len2 - 1.
 */

void _qadic_dense_mul_char_2(fmpz *rop, 
                             const fmpz *op1, long len1, const fmpz *op2, long len2, 
                             const fmpz *mod, const fmpz *invmod, long lenmod, long N)
{
    fmpz *t;

    t = _fmpz_vec_init(len1 + len2 - 1);

    /*_fmpz_mod_poly_mul(t, op1, len1, op2, len2, pN);*/
    _fmpz_poly_mul(t, op1, len1, op2, len2);
    _fmpz_vec_scalar_fdiv_r_2exp(t, t, len1 + len2 - 1, N);
    _qadic_dense_reduce_char_2(rop, t, len1 + len2 - 1, mod, invmod, lenmod, N);

    _fmpz_vec_clear(t, len1 + len2 - 1);
}

void _qadic_dense_mul(fmpz *rop, 
                const fmpz *op1, long len1, const fmpz *op2, long len2, 
                const fmpz *mod, const fmpz *invmod, long lenmod, const fmpz_t pN)
{
    fmpz *t;

    t = _fmpz_vec_init(len1 + len2 - 1);

    _fmpz_mod_poly_mul(t, op1, len1, op2, len2, pN);
    /* _fmpz_poly_mul(t, op1, len1, op2, len2); */
    _qadic_dense_reduce(rop, t, len1 + len2 - 1, mod, invmod, lenmod, pN);

    _fmpz_vec_clear(t, len1 + len2 - 1);
}

void qadic_dense_mul(qadic_dense_t x, const qadic_dense_t y, const qadic_dense_t z, 
                          const qadic_dense_ctx_t ctx)
{
    const long leny = y->length;
    const long lenz = z->length;
    const long lenx = leny + lenz - 1;
    const long N    = (&ctx->pctx)->N;
    const long d    = qadic_dense_ctx_degree(ctx);

    if (leny == 0 || lenz == 0 || y->val + z->val >= N)
    {
        qadic_dense_zero(x);
    }
    else if (!COEFF_IS_MPZ(*(&ctx->pctx)->p) && *(&ctx->pctx)->p == 2L)
    {
        fmpz *t, *mod, *invmod;

        x->val = y->val + z->val;

        if (x->val > 0)
        {
            mod = _fmpz_vec_init(d + 1);
            invmod = _fmpz_vec_init(d - 1);

            _fmpz_vec_scalar_fdiv_r_2exp(mod, ctx->mod->coeffs, d + 1, N - x->val);
            _fmpz_vec_scalar_fdiv_r_2exp(invmod, ctx->invmod->coeffs, d - 1, N - x->val);
        }
        else
        {
            mod = ctx->mod->coeffs;
            invmod = ctx->invmod->coeffs;
        }

        if (x == y || x == z)
        {
            t = _fmpz_vec_init(lenx);
        }
        else
        {
            padic_poly_fit_length(x, lenx);
            t = x->coeffs;
        }

        if (leny >= lenz)
            _qadic_dense_mul_char_2(t, y->coeffs, leny, z->coeffs, lenz,
                                    mod, invmod, d + 1, N - x->val);
        else
            _qadic_dense_mul_char_2(t, z->coeffs, lenz, y->coeffs, leny,
                                    mod, invmod, d + 1, N - x->val);

        if (x == y || x == z)
        {
            _fmpz_vec_clear(x->coeffs, x->alloc);
            x->coeffs = t;
            x->alloc  = lenx;
        }

        _padic_poly_set_length(x, FLINT_MIN(lenx, d));
        _padic_poly_normalise(x);

        if (x->val > 0)
        {
            _fmpz_vec_clear(mod, d + 1);
            _fmpz_vec_clear(invmod, d - 1);
        }
    }
    else
    {
        fmpz *t, *mod, *invmod;
        fmpz_t pN;
        int alloc;

        x->val = y->val + z->val;

        alloc = _padic_ctx_pow_ui(pN, N - x->val, &ctx->pctx);

        if (x->val > 0)
        {
            mod = _fmpz_vec_init(d + 1);
            invmod = _fmpz_vec_init(d - 1);

            _fmpz_vec_scalar_mod_fmpz(mod, ctx->mod->coeffs, d + 1, pN);
            _fmpz_vec_scalar_mod_fmpz(invmod, ctx->invmod->coeffs, d - 1, pN);
        }
        else
        {
            mod = ctx->mod->coeffs;
            invmod = ctx->invmod->coeffs;
        }

        if (x == y || x == z)
        {
            t = _fmpz_vec_init(lenx);
        }
        else
        {
            padic_poly_fit_length(x, lenx);
            t = x->coeffs;
        }

        if (leny >= lenz)
            _qadic_dense_mul(t, y->coeffs, leny, z->coeffs, lenz,
                             mod, invmod, d + 1, pN);
        else
            _qadic_dense_mul(t, z->coeffs, lenz, y->coeffs, leny,
                             mod, invmod, d + 1, pN);

        if (x == y || x == z)
        {
            _fmpz_vec_clear(x->coeffs, x->alloc);
            x->coeffs = t;
            x->alloc  = lenx;
        }

        _padic_poly_set_length(x, FLINT_MIN(lenx, d));
        _padic_poly_normalise(x);

        if (x->val > 0)
        {
            _fmpz_vec_clear(mod, d + 1);
            _fmpz_vec_clear(invmod, d - 1);
        }

        if (alloc)
            fmpz_clear(pN);
    }
}

