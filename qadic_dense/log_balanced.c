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

#include "fmpz_mod_poly.h"
#include "qadic_dense.h"

extern long _padic_log_bound(long v, long N, long p);

/*
    Assumes that P, T are vectors of length 2 d - 1.

    Assumes that 0 < len <= d.

    Assumes that 1 <= lo < hi.

    Does not support aliasing.
 */

static void 
_qadic_dense_log_bsplit_series(fmpz *P, fmpz_t B, fmpz *T, 
                         const fmpz *y, long len, long lo, long hi, 
                         const fmpz *mod, const fmpz *invmod, long lenmod)
{
    const long d = lenmod - 1;

    if (hi - lo == 1)
    {
        _fmpz_vec_set(P, y, len);
        _fmpz_vec_zero(P + len, 2*d - 1 - len);

        fmpz_set_si(B, lo);

        _fmpz_vec_set(T, P, 2*d - 1);
    }
    else if (hi - lo == 2)
    {
        _fmpz_poly_sqr(T, y, len);
        _fmpz_vec_zero(T + (2*len - 1), d - (2*len - 1));
        _qadic_dense_reduce_no_mod(P, T, 2*len - 1, mod, invmod, lenmod);

        fmpz_set_si(B, lo);
        fmpz_mul_si(B, B, lo + 1);

        _fmpz_vec_scalar_mul_si(T, y, len, lo + 1);
        _fmpz_vec_zero(T + len, d - len);
        _fmpz_vec_scalar_addmul_si(T, P, d, lo);

    }
    else
    {
        const long m = (lo + hi) / 2;

        fmpz *t, *RP, *RT, *W;
        fmpz_t RB;

        t  = _fmpz_vec_init(2*d - 1);
        RP = _fmpz_vec_init(2*d - 1);
        RT = _fmpz_vec_init(2*d - 1);
        W  = _fmpz_vec_init(2*d - 1);
        fmpz_init(RB);

        _qadic_dense_log_bsplit_series(P, B, T, y, len, lo, m, mod, invmod, lenmod);

        _qadic_dense_log_bsplit_series(RP, RB, RT, y, len, m, hi, mod, invmod, lenmod);

        _fmpz_poly_mul(t, RT, d, P, d);
        _qadic_dense_reduce_no_mod(W, t, 2*d - 1, mod, invmod, lenmod);
        _fmpz_vec_swap(RT, W, d);

        _fmpz_vec_scalar_mul_fmpz(T, T, d, RB);
        _fmpz_vec_scalar_addmul_fmpz(T, RT, d, B);

        _fmpz_poly_mul(t, P, d, RP, d);
        _qadic_dense_reduce_no_mod(W, t, 2*d - 1, mod, invmod, lenmod);
        _fmpz_vec_swap(P, W, d);

        fmpz_mul(B, B, RB);

        _fmpz_vec_clear(t,  2*d - 1);
        _fmpz_vec_clear(RP, 2*d - 1);
        _fmpz_vec_clear(RT, 2*d - 1);
        _fmpz_vec_clear(W,  2*d - 1);
        fmpz_clear(RB);
    }
}

/*
    Sets (z, d) to the sum 

        sum_{i=1}^{\infty} y^i / i mod p^N.

    The result may not be reduced modulo p^N, but it is 
    reduced modulo f(X) given by the data (a, j, lena).

    Supports aliasing between y and z.
 */

static void 
_qadic_dense_log_bsplit(fmpz *z, const fmpz *y, long v, long len, 
                  const fmpz *mod, const fmpz *invmod, long lenmod, 
                  const fmpz_t p, long N)
{
    const long d = lenmod - 1;

    fmpz *P, *T;
    fmpz_t B, C;
    long n;

    if (fmpz_fits_si(p))
    {
        n = _padic_log_bound(v, N, fmpz_get_si(p));
    }
    else
        n = (N - 1) / v;

    n = FLINT_MAX(n, 2);

    P = _fmpz_vec_init(2*d - 1);
    T = _fmpz_vec_init(2*d - 1);
    fmpz_init(B);
    fmpz_init(C);

    _qadic_dense_log_bsplit_series(P, B, T, y, len, 1, n, mod, invmod, lenmod);

    n = fmpz_remove(B, B, p);
    fmpz_pow_ui(C, p, n);
    _fmpz_vec_scalar_divexact_fmpz(T, T, d, C);

    _padic_inv(B, B, p, N);
    _fmpz_vec_scalar_mul_fmpz(z, T, d, B);

    _fmpz_vec_clear(P, 2*d - 1);
    _fmpz_vec_clear(T, 2*d - 1);
    fmpz_clear(B);
    fmpz_clear(C);
}

/*
    Computes 
    \begin{equation*}
    z = - \sum_{i = 1}^{\infty} \frac{y^i}{i} \pmod{p^N}.
    \end{equation*}

    Note that this can be used to compute the $p$-adic logarithm 
    via the equation 
    \begin{align*}
    \log(x) & = \sum_{i=1}^{\infty} (-1)^{i-1} \frac{(x-1)^i}{i} \\
            & = - \sum_{i=1}^{\infty} \frac{(1-x)^i}{i}.
    \end{align*}

    Assumes that $y = 1 - x$ is non-zero and that $v = \ord_p(y)$ 
    is at least $1$ when $p$ is odd and at least $2$ when $p = 2$ 
    so that the series converges.

    Assumes that $v < N$.

    Sets $(z, d)$.

    Ensures that the result is reduced modulo $p^N$.

    Does not support aliasing between $y$ and $z$.
 */

void _qadic_dense_log_balanced(fmpz *z, const fmpz *y, long len, 
                         const fmpz *mod, const fmpz *invmod, long lenmod, 
                         const fmpz_t p, long N, const fmpz_t pN)
{
    const long d = lenmod - 1;

    fmpz_t pv;
    fmpz *r, *s, *t, *u, *v;
    long i, w;

    r = _fmpz_vec_init(d);
    s = _fmpz_vec_init(2*d - 1);
    t = _fmpz_vec_init(d);
    u = _fmpz_vec_init(d);
    v = _fmpz_vec_init(2*d - 1);
    fmpz_init(pv);

    fmpz_set(pv, p);
    _fmpz_vec_scalar_mod_fmpz(t, y, len, pN);
    _fmpz_vec_zero(z, d);
    w = 1;

    while (!_fmpz_vec_is_zero(t, d))
    {
        fmpz_mul(pv, pv, pv);

        for (i = 0; i < d; i++)
        {
            fmpz_fdiv_qr(t + i, r + i, t + i, pv);
        }

        if (!_fmpz_vec_is_zero(t, d))
        {
            _fmpz_vec_scalar_mul_fmpz(t, t, d, pv);

            /* Temporarily set r = 1 - r to compute u = (1-r)^{-1} */
            fmpz_sub_ui(r + 0, r + 0, 1);
            _fmpz_vec_neg(r, r, d);
            _qadic_dense_inv(u, r, d, mod, invmod, lenmod, p, N);
            _fmpz_vec_neg(r, r, d);
            fmpz_add_ui(r + 0, r + 0, 1);

            _fmpz_poly_mul(v, t, d, u, d);
            _qadic_dense_reduce_no_mod(s, v, 2 * d - 1, mod, invmod, lenmod);
            _fmpz_vec_scalar_mod_fmpz(t, s, d, pN);
        }

        if (!_fmpz_vec_is_zero(r, d))
        {
            _qadic_dense_log_bsplit(r, r, w, d, mod, invmod, lenmod, p, N);
            _fmpz_vec_sub(z, z, r, d);
            _fmpz_vec_scalar_mod_fmpz(z, z, d, pN);
        }

        w *= 2;
    }

    _fmpz_vec_clear(r, d);
    _fmpz_vec_clear(s, 2*d - 1);
    _fmpz_vec_clear(t, d);
    _fmpz_vec_clear(u, d);
    _fmpz_vec_clear(v, 2*d - 1);
    fmpz_clear(pv);
}

int qadic_dense_log_balanced(qadic_dense_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx)
{
    const fmpz *p  = (&ctx->pctx)->p;
    const long d   = qadic_dense_ctx_degree(ctx);
    const long N   = (&ctx->pctx)->N;
    const long len = op->length;

    if (op->val < 0)
    {
        return 0;
    }
    else
    {
        fmpz *x;
        fmpz_t pN;
        int alloc, ans;

        x = _fmpz_vec_init(len + 1);
        alloc = _padic_ctx_pow_ui(pN, N, &ctx->pctx);

        /* Set x := (1 - op) mod p^N */
        fmpz_pow_ui(x + len, p, op->val);
        _fmpz_vec_scalar_mul_fmpz(x, op->coeffs, len, x + len);
        fmpz_sub_ui(x + 0, x + 0, 1);
        _fmpz_vec_neg(x, x, len);
        _fmpz_vec_scalar_mod_fmpz(x, x, len, pN);

        if (_fmpz_vec_is_zero(x, len))
        {
            padic_poly_zero(rop);
            ans = 1;
        }
        else
        {
            const long v = _fmpz_vec_ord_p(x, len, p);

            if (v >= 2 || (*p != 2L && v >= 1))
            {
                if (v >= N)
                {
                    padic_poly_zero(rop);
                }
                else
                {
                    padic_poly_fit_length(rop, d);

                    _qadic_dense_log_balanced(rop->coeffs, x, len,
                                              ctx->mod->coeffs, ctx->invmod->coeffs, d + 1,
                                              p, N, pN);
                    rop->val = 0;

                    _padic_poly_set_length(rop, d);
                    _padic_poly_normalise(rop);
                    padic_poly_canonicalise(rop, p);
                }
                ans = 1;
            }
            else
            {
                ans = 0;
            }
        }

        _fmpz_vec_clear(x, len + 1);
        if (alloc)
            fmpz_clear(pN);
        return ans;
    }
}

