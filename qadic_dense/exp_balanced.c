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

extern long _padic_exp_bound(long v, long N, const fmpz_t p);

static void 
_qadic_dense_exp_bsplit_series(fmpz *P, fmpz_t Q, fmpz *T,
                               const fmpz *x, long len, long lo, long hi,
                               const fmpz *mod, const fmpz *invmod, long lenmod)
{
    const long d = lenmod - 1;

    if (hi - lo == 1)
    {
        _fmpz_vec_set(P, x, len);
        _fmpz_vec_zero(P + len, 2*d - 1 - len);

        fmpz_set_si(Q, lo);

        _fmpz_vec_set(T, P, 2*d - 1);
    }
    else if (hi - lo == 2)
    {
        _fmpz_poly_sqr(T, x, len);
        _fmpz_vec_zero(T + (2*len - 1), d - (2*len - 1));
        _qadic_dense_reduce_no_mod(P, T, 2*len - 1, mod, invmod, lenmod);

        fmpz_set_si(Q, lo);
        fmpz_mul_si(Q, Q, lo + 1);

        _fmpz_vec_scalar_mul_si(T, x, len, lo + 1);
        _fmpz_vec_zero(T + len, d - len);
        _fmpz_vec_add(T, T, P, d);
    }
    else
    {
        const long m = (lo + hi) / 2;

        fmpz *t, *PR, *TR, *W;
        fmpz_t QR;

        t  = _fmpz_vec_init(2*d - 1);
        PR = _fmpz_vec_init(2*d - 1);
        TR = _fmpz_vec_init(2*d - 1);
        W  = _fmpz_vec_init(2*d - 1);
        fmpz_init(QR);

        _qadic_dense_exp_bsplit_series(P, Q, T, x, len, lo, m, mod, invmod, lenmod);

        _qadic_dense_exp_bsplit_series(PR, QR, TR, x, len, m, hi, mod, invmod, lenmod);

        _fmpz_poly_mul(t, TR, d, P, d);
        _qadic_dense_reduce_no_mod(W, t, 2*d - 1, mod, invmod, lenmod);

        _fmpz_vec_scalar_mul_fmpz(T, T, d, QR);
        _fmpz_vec_add(T, T, W, d);

        _fmpz_poly_mul(t, P, d, PR, d);
        _qadic_dense_reduce_no_mod(W, t, 2*d - 1, mod, invmod, lenmod);
        _fmpz_vec_swap(P, W, d);

        fmpz_mul(Q, Q, QR);

        _fmpz_vec_clear(t,  2*d - 1);
        _fmpz_vec_clear(PR, 2*d - 1);
        _fmpz_vec_clear(TR, 2*d - 1);
        _fmpz_vec_clear(W,  2*d - 1);
        fmpz_clear(QR);
    }
}

static void 
_qadic_dense_exp_bsplit(fmpz *y, const fmpz *x, long v, long len, 
                  const fmpz *mod, const fmpz* invmod, long lenmod, 
                  const fmpz_t p, long N)
{
    const long d = lenmod - 1;
    const long n = _padic_exp_bound(v, N, p);

    if (n == 1)
    {
        fmpz_one(y + 0);
        _fmpz_vec_zero(y + 1, d - 1);
    }
    else
    {
        fmpz *P, *T;
        fmpz_t Q, R;
        long f;

        P = _fmpz_vec_init(2*d - 1);
        T = _fmpz_vec_init(2*d - 1);
        fmpz_init(Q);
        fmpz_init(R);

        _qadic_dense_exp_bsplit_series(P, Q, T, x, len, 1, n, mod, invmod, lenmod);

        fmpz_add(T + 0, T + 0, Q);  /* (T,Q) := (T,Q) + 1 */

        /* Note exp(x) is a unit so val(T) == val(Q). */
        f = fmpz_remove(Q, Q, p);
        fmpz_pow_ui(R, p, f);
        _fmpz_vec_scalar_divexact_fmpz(T, T, d, R);

        _padic_inv(Q, Q, p, N);
        _fmpz_vec_scalar_mul_fmpz(y, T, d, Q);

        _fmpz_vec_clear(P, 2*d - 1);
        _fmpz_vec_clear(T, 2*d - 1);
        fmpz_clear(Q);
        fmpz_clear(R);
    }
}

void _qadic_dense_exp_balanced(fmpz *rop, const fmpz *x, long v, long len, 
                         const fmpz *mod, const fmpz *invmod, long lenmod, 
                         const fmpz_t p, long N, const fmpz_t pN)
{
    const long d = lenmod - 1;

    fmpz_t pw;
    fmpz *r, *s, *t;
    long i, w;

    r = _fmpz_vec_init(2*d - 1);
    s = _fmpz_vec_init(2*d - 1);
    t = _fmpz_vec_init(d);
    fmpz_init(pw);

    fmpz_pow_ui(pw, p, v);
    _fmpz_vec_scalar_mul_fmpz(t, x, len, pw);
    _fmpz_vec_scalar_mod_fmpz(t, t, len, pN);
    _fmpz_vec_zero(t + len, d - len);

    fmpz_set(pw, p);
    fmpz_one(rop + 0);
    _fmpz_vec_zero(rop + 1, d - 1);
    w = 1;

    while (!_fmpz_vec_is_zero(t, d))
    {
        fmpz_mul(pw, pw, pw);

        for (i = 0; i < d; i++)
        {
            fmpz_fdiv_r(r + i, t + i, pw);
            fmpz_sub(t + i, t + i, r + i);
        }

        if (!_fmpz_vec_is_zero(r, d))
        {
            _qadic_dense_exp_bsplit(r, r, w, d, mod, invmod, lenmod, p, N);
            _fmpz_poly_mul(s, rop, d, r, d);
            _qadic_dense_reduce_no_mod(r, s, 2*d - 1, mod, invmod, lenmod);
            _fmpz_vec_scalar_mod_fmpz(rop, r, d, pN);
        }

        w *= 2;
    }

    _fmpz_vec_clear(r, 2*d - 1);
    _fmpz_vec_clear(s, 2*d - 1);
    _fmpz_vec_clear(t, d);
    fmpz_clear(pw);
}

int qadic_dense_exp_balanced(qadic_dense_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx)
{
    const long N  = (&ctx->pctx)->N;
    const long v  = op->val;
    const fmpz *p = (&ctx->pctx)->p;

    if (padic_poly_is_zero(op))
    {
        padic_poly_one(rop, &ctx->pctx);
        return 1;
    }

    if ((*p == 2L && v <= 1) || (v <= 0))
    {
        return 0;
    }
    else
    {
        if (v < N)
        {
            const long d = qadic_dense_ctx_degree(ctx);

            fmpz_t pN;
            int alloc;

            alloc = _padic_ctx_pow_ui(pN, N, &ctx->pctx);

            padic_poly_fit_length(rop, d);

            _qadic_dense_exp_balanced(rop->coeffs, op->coeffs, v, op->length, 
                                      ctx->mod->coeffs, ctx->invmod->coeffs, d + 1, p, N, pN);
            rop->val = 0;

            _padic_poly_set_length(rop, d);
            _padic_poly_normalise(rop);

            if (alloc)
                fmpz_clear(pN);
        }
        else
        {
            padic_poly_one(rop, &ctx->pctx);
        }
        return 1;
    }
}

