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
#include "ulong_extras.h"
#include "qadic_dense.h"

/*
    Assumes that \code{len1} and \code{len2} are positive but at 
    most~$d$, and also that \code{len1} is at least $6$.

    The latter assumption guarantees that $\ceil{n/B} \geq 2$, 
    i.e.\ $n \geq 2B$ so $n \geq 2 \ceil{\sqrt{n}}$.
 */

static void 
_qadic_dense_compose_mod_rectangular(fmpz *rop, 
                           const fmpz *op1, long len1, 
                           const fmpz *op2, long len2, 
                           const fmpz *mod, const fmpz *invmod, long lenmod, 
                           const fmpz_t p)
{
    const long d = lenmod - 1;

    if (len2 == 1)
    {
        _fmpz_mod_poly_evaluate_fmpz(rop, op1, len1, op2, p);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else
    {
        const long B = n_sqrt(len1);
        long i, k;
        fmpz *pows, *s, *t;

        pows = _fmpz_vec_init((B + 2) * d);
        s    = _fmpz_vec_init(2 * d - 1);
        t    = _fmpz_vec_init(2 * d - 1);

        fmpz_one(pows + 0 * d + 0);
        _fmpz_vec_set(pows + 1 * d, op2, len2);
        for (i = 2; i <= B; i++)
        {
            _fmpz_poly_mul(t, pows + (i - 1) * d, d, op2, len2);
            _qadic_dense_reduce_no_mod(s, t, d + len2 - 1, mod, invmod, lenmod);
            _fmpz_vec_scalar_mod_fmpz(pows + i * d, s, d, p);
        }

        _fmpz_vec_zero(rop, d);

        for (i = (len1 + B - 1) / B - 1; i >= 0; i--)
        {
            _fmpz_poly_mul(s, rop, d, pows + B * d, d);
            _qadic_dense_reduce(t, s, 2 * d - 1, mod, invmod, lenmod, p);

            _fmpz_vec_set(rop, t, d);
            fmpz_add(rop + 0, rop + 0, op1 + i*B);
            for (k = FLINT_MIN(B, len1 - i*B) - 1; k > 0; k--)
            {
                _fmpz_vec_scalar_addmul_fmpz(rop, pows + k * d, d, op1 + (i*B + k));
            }

            _fmpz_vec_scalar_mod_fmpz(rop, rop, d, p);
        }

        _fmpz_vec_clear(pows, (B + 2) * d);
        _fmpz_vec_clear(s, 2 * d - 1);
        _fmpz_vec_clear(t, 2 * d - 1);
    }
}

static void 
_qadic_dense_compose_mod_horner(fmpz *rop, 
                           const fmpz *op1, long len1, 
                           const fmpz *op2, long len2, 
                           const fmpz *mod, const fmpz *invmod, long lenmod, 
                           const fmpz_t p)
{
    const long d = lenmod - 1;

    if (len1 == 1)
    {
        fmpz_set(rop, op1);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else if (len2 == 1)
    {
        _fmpz_mod_poly_evaluate_fmpz(rop, op1, len1, op2, p);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else
    {
        long i;
        fmpz *s, *t;

        t = _fmpz_vec_init(2*d - 1);
        s = _fmpz_vec_init(2*d - 1);

        _fmpz_vec_zero(rop, d);

        for (i = len1 - 1; i >= 0; i--)
        {
            _fmpz_poly_mul(s, rop, d, op2, len2);
            _qadic_dense_reduce_no_mod(t, s, d + len2 - 1, mod, invmod, lenmod);
            _fmpz_poly_add(rop, t, d, op1 + i, 1);
            _fmpz_vec_scalar_mod_fmpz(rop, rop, d, p);
        }

        _fmpz_vec_clear(s, 2*d - 1);
        _fmpz_vec_clear(t, 2*d - 1);
    }
}

/* 
    Computes the composition $f(g(X))$ modulo the sparse polynomial 
    given by the data \code{(a, j, lena)}, which is assumed to be 
    of degree~$d \geq 2$.

    Sets the vector \code{(rop, d)}.

    Assumes that \code{len1} and \code{len2} are positive but at 
    most~$d$.

    Does not support aliasing.
 */

void 
_qadic_dense_compose_mod(fmpz *rop, 
                           const fmpz *op1, long len1, 
                           const fmpz *op2, long len2, 
                           const fmpz *mod, const fmpz *invmod, long lenmod, 
                           const fmpz_t p)
{
    if (len1 < 6)
    {
        _qadic_dense_compose_mod_horner(rop, op1, len1, op2, len2, mod, invmod, lenmod, p);
    }
    else
    {
        _qadic_dense_compose_mod_rectangular(rop, op1, len1, op2, len2, mod, invmod, lenmod, p);
    }
}

void _qadic_dense_frobenius_a(fmpz *rop, long exp, 
                    const fmpz *mod, const fmpz *invmod, long lenmod, 
                    const fmpz_t p, long N)
{
    const long d = lenmod - 1;

    long *e, i, n;
    fmpz *pow, *f1, *f2, *inv, *s, *t;

    n = FLINT_CLOG2(N) + 1;

    e = flint_malloc(n * sizeof(long));
    for (e[i = 0] = N; e[i] > 1; i++)
        e[i + 1] = (e[i] + 1) / 2;

    pow = _fmpz_vec_init(n);
    f1  = _fmpz_vec_init(d + 1);
    f2  = _fmpz_vec_init(d);
    inv = _fmpz_vec_init(2*d - 1);
    s   = _fmpz_vec_init(2*d - 1);
    t   = _fmpz_vec_init(2*d - 1);

    /* Compute powers of p */
    {
        fmpz_one(t);
        fmpz_set(pow + i, p);
    }
    for (i--; i >= 1; i--)
    {
        if (e[i] & 1L)
        {
            fmpz_mul(pow + i, t, pow + (i + 1));
            fmpz_mul(t, t, t);
        }
        else
        {
            fmpz_mul(t, t, pow + (i + 1));
            fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
        }
    }
    {
        if (e[i] & 1L)
            fmpz_mul(pow + i, t, pow + (i + 1));
        else
            fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
    }

    /* Dense representation of f and f' */
    {
        long k;

        _fmpz_vec_set(f1, mod, lenmod);
        for (k = 1; k < lenmod; k++)
            fmpz_mul_ui(f2 + (k - 1), mod + k, k);
    }

    /* Run Newton iteration */
    i = n - 1;
    {
        fmpz op[2] = {0L, 1L};

        fmpz_pow_ui(t, p, exp);
        _qadic_dense_pow(rop, op, 2, t, mod, invmod, lenmod, pow + i);
        _qadic_dense_compose_mod(t, f2, d, rop, d, mod, invmod, lenmod, pow + i);
        _qadic_dense_inv(inv, t, d, mod, invmod, lenmod, p, 1);
    }
    for (i--; i >= 0; i--)
    {
        _qadic_dense_compose_mod(s, f1, d + 1, rop, d, mod, invmod, lenmod, pow + i);
        _fmpz_mod_poly_mul(t, s, d, inv, d, pow + i);
        _qadic_dense_reduce(s, t, 2*d - 1, mod, invmod, lenmod, pow + i);
        _fmpz_mod_poly_sub(rop, rop, d, s, d, pow + i);

        if (i > 0)
        {
            _qadic_dense_compose_mod(s, f2, d, rop, d, mod, invmod, lenmod, pow + i);
            _fmpz_mod_poly_mul(t, inv, d, s, d, pow + i);
            _qadic_dense_reduce(s, t, 2*d - 1, mod, invmod, lenmod, pow + i);
            fmpz_sub_ui(s, s, 2);
            if (fmpz_sgn(s) < 0)
                fmpz_add(s, s, pow + i);
            _fmpz_mod_poly_neg(s, s, d, pow + i);
            _fmpz_mod_poly_mul(t, inv, d, s, d, pow + i);
            _qadic_dense_reduce(s, t, 2*d - 1, mod, invmod, lenmod, pow + i);

            /* SWAP(inv, s).  Requires the arrays to be of the same size. */
            {
                fmpz *__t;

                __t = inv;
                inv = s;
                s   = __t;
            }
        }
    }

    _fmpz_vec_clear(pow, n);
    _fmpz_vec_clear(f1, d + 1);
    _fmpz_vec_clear(f2, d);
    _fmpz_vec_clear(inv, 2*d - 1);
    _fmpz_vec_clear(s, 2*d - 1);
    _fmpz_vec_clear(t, 2*d - 1);
    flint_free(e);
}

/*
    Sets (rop, 2d-1) to the image of (op, len) under the Frobenius operator 
    raised to the e-th power.
 */

void _qadic_dense_frobenius(fmpz *rop, const fmpz *op, long len, long e, 
                  const fmpz *mod, const fmpz *invmod, long lenmod, 
                  const fmpz_t p, long N)
{
    const long d = lenmod - 1;

    if (len == 1)  /* op is in Zp, not just Zq */
    {
        _fmpz_vec_set(rop, op, len);
        _fmpz_vec_zero(rop + len, (2*d - 1)  - len);
    }
    else if (N == 1)
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_pow_ui(t, p, e);
        _qadic_dense_pow(rop, op, len, t, mod, invmod, lenmod, p);
        fmpz_clear(t);
    }
    else
    {
        fmpz *t;
        fmpz_t pow;

        t = _fmpz_vec_init(2*d - 1);
        fmpz_init(pow);
        fmpz_pow_ui(pow, p, N);

        _qadic_dense_frobenius_a(t, e, mod, invmod, lenmod, p, N);

        _qadic_dense_compose_mod(rop, op, len, t, d, mod, invmod, lenmod, pow);
        _fmpz_vec_zero(rop + d, d - 1);

        _fmpz_vec_clear(t, 2*d - 1);
        fmpz_clear(pow);
    }
}

void qadic_dense_frobenius(qadic_dense_t rop, const qadic_dense_t op, long e, const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;
    const long d = qadic_dense_ctx_degree(ctx);

    e = e % d;
    if (e < 0)
        e += d;

    if (qadic_dense_is_zero(op) || op->val >= N)
    {
        qadic_dense_zero(rop);
    }
    else if (e == 0)
    {
        padic_poly_set(rop, op);
        padic_poly_reduce(rop, &ctx->pctx);
    }
    else
    {
        fmpz *t;
        fmpz *mod, *invmod;

        if (op->val >= 0)
        {
            int alloc;
            fmpz_t pow;

            alloc = _padic_ctx_pow_ui(pow, N - op->val, &ctx->pctx);

            mod = _fmpz_vec_init(d + 1);
            _fmpz_vec_scalar_mod_fmpz(mod, ctx->mod->coeffs, d + 1, pow);
            invmod = _fmpz_vec_init(d + 1);
            _fmpz_vec_scalar_mod_fmpz(invmod, ctx->invmod->coeffs, d + 1, pow);

            if (alloc)
                fmpz_clear(pow);
        }
        else
        {
            mod = ctx->mod->coeffs;
            invmod = _fmpz_vec_init(d + 1);
            _qadic_dense_ctx_init_inv(invmod, mod, &ctx->pctx, d, N - op->val);
        }

        if (rop == op)
        {
            t = _fmpz_vec_init(2 * d - 1);
        }
        else
        {
            padic_poly_fit_length(rop, 2 * d - 1);
            t = rop->coeffs;
        }

        _qadic_dense_frobenius(t, op->coeffs, op->length, e,
                               mod, invmod, d + 1,
                               (&ctx->pctx)->p, N - op->val);

        if (op->val >= 0)
        {
            flint_free(mod);
            flint_free(invmod);
        }
        else
        {
            flint_free(invmod);
        }

        if (rop == op)
        {
            _fmpz_vec_clear(rop->coeffs, rop->alloc);
            rop->coeffs = t;
            rop->alloc  = 2 * d - 1;
            rop->length = d;
        }
        else
        {
            rop->val = op->val;
            _padic_poly_set_length(rop, d);
        }
        _padic_poly_normalise(rop);
    }
}

