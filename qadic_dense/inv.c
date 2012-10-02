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

#include "fmpz_mod_poly.h"
#include "qadic_dense.h"

void _qadic_dense_inv(fmpz *rop, const fmpz *op, long len, 
                const fmpz *mod, const fmpz *invmod, long lenmod, 
                const fmpz_t p, long N)
{
    const long d = lenmod - 1;

    if (len == 1)
    {
        _padic_inv(rop, op, p, N);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else if (N == 1)
    {
        fmpz *P = _fmpz_vec_init(lenmod);

        _fmpz_vec_scalar_mod_fmpz(P, mod, lenmod, p);

        _fmpz_mod_poly_invmod(rop, op, len, P, lenmod, p);

        _fmpz_vec_clear(P, lenmod);
    }
    else  /* d, N >= 2 */
    {
        long *e, i, n;
        fmpz *pow, *redmod, *redinvmod, *u;
        fmpz *s, *t;

        n = FLINT_CLOG2(N) + 1;

        /* Compute sequence of exponents */
        e = flint_malloc(n * sizeof(long));
        for (e[i = 0] = N; e[i] > 1; i++)
            e[i + 1] = (e[i] + 1) / 2;

        pow       = _fmpz_vec_init(n);
        redmod    = _fmpz_vec_init(lenmod * n);
        redinvmod = _fmpz_vec_init(lenmod * n);
        u = _fmpz_vec_init(len * n);
        s = _fmpz_vec_init(2 * d - 1);
        t = _fmpz_vec_init(2 * d - 1);

        /* Compute powers of p and reduced moduli */
        {
            fmpz_one(t);
            fmpz_set(pow + i, p);
            _fmpz_vec_scalar_mod_fmpz(redmod + i * lenmod, mod, lenmod, pow + i);
            _fmpz_vec_scalar_mod_fmpz(redinvmod + i * lenmod, invmod, lenmod, pow + i);
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
            _fmpz_vec_scalar_mod_fmpz(redmod + i * lenmod, mod, lenmod, pow + i);
            _fmpz_vec_scalar_mod_fmpz(redinvmod + i * lenmod, invmod, lenmod, pow + i);
        }
        {
            if (e[i] & 1L)
                fmpz_mul(pow + i, t, pow + (i + 1));
            else
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
            _fmpz_vec_scalar_mod_fmpz(redmod + i * lenmod, mod, lenmod, pow + i);
            _fmpz_vec_scalar_mod_fmpz(redinvmod + i * lenmod, invmod, lenmod, pow + i);
        }

        /* Compute reduced units */
        {
            _fmpz_vec_scalar_mod_fmpz(u + 0 * len, op, len, pow + 0);
        }
        for (i = 1; i < n; i++)
        {
            _fmpz_vec_scalar_mod_fmpz(u + i * len, u + (i - 1) * len, len, pow + i);
        }

        /* Run Newton iteration */
        i = n - 1;
        {
            _fmpz_mod_poly_invmod(rop, u + i * len, len, redmod + i * lenmod, lenmod, pow + i);
        }
        for (i--; i >= 0; i--)  /* z' := 2 z - a z^2 */
        {
            _fmpz_mod_poly_sqr(s, rop, d, pow + i);
            _qadic_dense_reduce(t, s, 2 * d - 1, redmod + i * lenmod,
                                redinvmod + i * lenmod, lenmod, pow + i);

            _fmpz_mod_poly_mul(s, t, d, u + i * len, len, pow + i);
            _qadic_dense_reduce(t, s, d + len - 1, redmod + i * lenmod,
                                redinvmod + i * lenmod, lenmod, pow + i);

            _fmpz_vec_scalar_mul_2exp(rop, rop, d, 1);
            _fmpz_mod_poly_sub(rop, rop, d, t, d, pow + i);
            /* _fmpz_poly_sub(rop, rop, d, t, d);
               _fmpz_vec_scalar_mod_fmpz(rop, rop, d, pow + i); */
        }

        _fmpz_vec_clear(pow, n);
        _fmpz_vec_clear(redmod, lenmod * n);
        _fmpz_vec_clear(redinvmod, lenmod * n);
        _fmpz_vec_clear(u, len * n);
        _fmpz_vec_clear(s, 2 * d - 1);
        _fmpz_vec_clear(t, 2 * d - 1);
        flint_free(e);
    }
}

void qadic_dense_inv(qadic_dense_t x, const qadic_dense_t y, const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (qadic_dense_is_zero(y))
    {
        printf("Exception (qadic_dense_inv).  Zero is not invertible.\n");
        abort();
    }

    /*
        If y = u p^v has negative valuation with N <= -v then the 
        exact inverse of y is zero when reduced modulo $p^N$
     */
    if (N + y->val <= 0)
    {
        qadic_dense_zero(x);
    }
    else
    {
        const long d = qadic_dense_ctx_degree(ctx);
        fmpz *t;

        if (x == y)
        {
            t = _fmpz_vec_init(d);
        }
        else
        {
            padic_poly_fit_length(x, d);
            t = x->coeffs;
        }

        _qadic_dense_inv(t, y->coeffs, y->length,
                         ctx->mod->coeffs, ctx->invmod->coeffs,
                         d + 1, (&ctx->pctx)->p, N + y->val);
        x->val = - y->val;

        if (x == y)
        {
            _fmpz_vec_clear(x->coeffs, x->alloc);
            x->coeffs = t;
            x->alloc  = d;
            x->length = d;
        }
        else
        {
            _padic_poly_set_length(x, d);
        }
        _padic_poly_normalise(x);
    }
}

