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

void _qadic_dense_pow_char_2(fmpz *rop, const fmpz *op, long len, const fmpz_t e, 
                   const fmpz *mod, const fmpz *invmod, long lenmod, long N)
{
    const long d = lenmod - 1;

    if (fmpz_is_zero(e))
    {
        fmpz_one(rop);
        _fmpz_vec_zero(rop + 1, 2 * d - 1 - 1);
    }
    else if (fmpz_is_one(e))
    {
        _fmpz_vec_set(rop, op, len);
        _fmpz_vec_zero(rop + len, 2 * d - 1 - len);
    }
    else
    {
        ulong bit;
        long l;
        fmpz *u = _fmpz_vec_init(2 * d - 1);
        fmpz *v = _fmpz_vec_init(2 * d - 1);
        fmpz *R, *S, *T;

        _fmpz_vec_zero(rop, 2 * d - 1);

        /*
           Set bits to the bitmask with a 1 one place lower than the msb of e
         */

        bit = fmpz_bits(e) - 2;

        /*
           Trial run without any polynomial arithmetic to determine the parity 
           of the number of swaps;  then set R and S accordingly
         */
        
        {
            unsigned int swaps = 0U;
            ulong bit2 = bit;
            if (fmpz_tstbit(e, bit2))
                swaps = ~swaps;
            while (bit2--)
                if (!fmpz_tstbit(e, bit2))
                    swaps = ~swaps;
            
            if (swaps == 0U)
            {
                R = rop;
                S = v;
            }
            else
            {
                R = v;
                S = rop;
            }
        }
        
        /*
           We unroll the first step of the loop, referring to {op, len}
         */

        _fmpz_poly_sqr(u, op, len);
        _fmpz_vec_scalar_fdiv_r_2exp(u, u, 2 * len - 1, N);
        _qadic_dense_reduce_char_2(R, u, 2 * len - 1, mod, invmod, lenmod, N);
        l = FLINT_MIN(2 * len - 1, d);
        FMPZ_VEC_NORM(R, l);

        if (fmpz_tstbit(e, bit))
        {
            if (l >= len)
                _fmpz_poly_mul(u, R, l, op, len);
            else
                _fmpz_poly_mul(u, op, len, R, l);
            _fmpz_vec_scalar_fdiv_r_2exp(u, u, l + len - 1, N);
            _qadic_dense_reduce_char_2(S, u, l + len - 1, mod, invmod, lenmod, N);
            l = FLINT_MIN(l + len - 1, d);
            FMPZ_VEC_NORM(S, l);
            T = R;
            R = S;
            S = T;
        }

        while (bit--)
        {
            if (fmpz_tstbit(e, bit))
            {
                _fmpz_poly_sqr(u, R, l);
                _fmpz_vec_scalar_fdiv_r_2exp(u, u, 2 * l - 1, N);
                _qadic_dense_reduce_char_2(S, u, 2 * l - 1, mod, invmod, lenmod, N);
                l = FLINT_MIN(2* l - 1, d);
                FMPZ_VEC_NORM(S, l);
                if (l >= len)
                    _fmpz_poly_mul(u, S, l, op, len);
                else
                    _fmpz_poly_mul(u, op, len, S, l);
                _fmpz_vec_scalar_fdiv_r_2exp(u, u, l + len - 1, N);
                _qadic_dense_reduce_char_2(R, u, l + len - 1, mod, invmod, lenmod, N);
                l = FLINT_MIN(l + len - 1, d);
                FMPZ_VEC_NORM(R, l);
            }
            else
            {
                 _fmpz_poly_sqr(u, R, l);
                _fmpz_vec_scalar_fdiv_r_2exp(u, u, 2 * l - 1, N);
                _qadic_dense_reduce_char_2(S, u, 2 * l - 1, mod, invmod, lenmod, N);
                l = FLINT_MIN(2 * l - 1, d);
                FMPZ_VEC_NORM(S, l);
                T = R;
                R = S;
                S = T;
            }
        }

        _fmpz_vec_clear(u, 2 * d - 1);
        _fmpz_vec_clear(v, 2 * d - 1);
    }
}

void _qadic_dense_pow(fmpz *rop, const fmpz *op, long len, const fmpz_t e, 
                   const fmpz *mod, const fmpz *invmod, long lenmod, const fmpz_t p)
{
    const long d = lenmod - 1;

    if (fmpz_is_zero(e))
    {
        fmpz_one(rop);
        _fmpz_vec_zero(rop + 1, 2 * d - 1 - 1);
    }
    else if (fmpz_is_one(e))
    {
        _fmpz_vec_set(rop, op, len);
        _fmpz_vec_zero(rop + len, 2 * d - 1 - len);
    }
    else
    {
        ulong bit;
        long l;
        fmpz *u = _fmpz_vec_init(2 * d - 1);
        fmpz *v = _fmpz_vec_init(2 * d - 1);
        fmpz *R, *S, *T;

        _fmpz_vec_zero(rop, 2 * d - 1);

        /*
           Set bits to the bitmask with a 1 one place lower than the msb of e
         */

        bit = fmpz_bits(e) - 2;

        /*
           Trial run without any polynomial arithmetic to determine the parity 
           of the number of swaps;  then set R and S accordingly
         */
        
        {
            unsigned int swaps = 0U;
            ulong bit2 = bit;
            if (fmpz_tstbit(e, bit2))
                swaps = ~swaps;
            while (bit2--)
                if (!fmpz_tstbit(e, bit2))
                    swaps = ~swaps;
            
            if (swaps == 0U)
            {
                R = rop;
                S = v;
            }
            else
            {
                R = v;
                S = rop;
            }
        }
        
        /*
           We unroll the first step of the loop, referring to {op, len}
         */

        _fmpz_mod_poly_sqr(u, op, len, p);
        /* _fmpz_poly_sqr(u, op, len); */
        _qadic_dense_reduce(R, u, 2 * len - 1, mod, invmod, lenmod, p);
        l = FLINT_MIN(2 * len - 1, d);
        FMPZ_VEC_NORM(R, l);

        if (fmpz_tstbit(e, bit))
        {
            if (l >= len)
                _fmpz_mod_poly_mul(u, R, l, op, len, p);
            else
                _fmpz_mod_poly_mul(u, op, len, R, l, p);
            /* _fmpz_poly_mul(u, R, d, op, len); */
            _qadic_dense_reduce(S, u, l + len - 1, mod, invmod, lenmod, p);
            l = FLINT_MIN(l + len - 1, d);
            FMPZ_VEC_NORM(S, l);
            T = R;
            R = S;
            S = T;
        }

        while (bit--)
        {
            if (fmpz_tstbit(e, bit))
            {
                _fmpz_mod_poly_sqr(u, R, l, p);
                /* _fmpz_poly_sqr(u, R, d); */
                _qadic_dense_reduce(S, u, 2 * l - 1, mod, invmod, lenmod, p);
                l = FLINT_MIN(2* l - 1, d);
                FMPZ_VEC_NORM(S, l);
                if (l >= len)
                    _fmpz_mod_poly_mul(u, S, l, op, len, p);
                else
                    _fmpz_mod_poly_mul(u, op, len, S, l, p);
                /* _fmpz_poly_mul(u, S, d, op, len); */
                _qadic_dense_reduce(R, u, l + len - 1, mod, invmod, lenmod, p);
                l = FLINT_MIN(l + len - 1, d);
                FMPZ_VEC_NORM(R, l);
            }
            else
            {
                 _fmpz_mod_poly_sqr(u, R, l, p);
                 /* _fmpz_poly_sqr(u, R, d); */
                _qadic_dense_reduce(S, u, 2 * l - 1, mod, invmod, lenmod, p);
                l = FLINT_MIN(2 * l - 1, d);
                FMPZ_VEC_NORM(S, l);
                T = R;
                R = S;
                S = T;
            }
        }

        _fmpz_vec_clear(u, 2 * d - 1);
        _fmpz_vec_clear(v, 2 * d - 1);
    }
}

void qadic_dense_pow(qadic_dense_t x, const qadic_dense_t y, const fmpz_t e, const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (fmpz_sgn(e) < 0)
    {
        printf("Exception (qadic_dense_pow).  e < 0.\n");
        abort();
    }

    if (fmpz_is_zero(e))
    {
        qadic_dense_one(x, ctx);
    }
    else if (qadic_dense_is_zero(y))
    {
        qadic_dense_zero(x);
    }
    else
    {
        fmpz_t val;  /* N - e * val(y) */

        fmpz_init_set(val, e);
        fmpz_mul_si(val, val, y->val);

        if (fmpz_cmp_si(val, N) >= 0)
        {
            qadic_dense_zero(x);
        }
        else if (fmpz_is_one(e))
        {
            qadic_dense_set(x, y);
            /*qadic_dense_reduce(x, ctx);*/
        }
        else if (!COEFF_IS_MPZ(*(&ctx->pctx)->p) && *(&ctx->pctx)->p == 2L)
        {
            const long d = qadic_dense_ctx_degree(ctx);
            fmpz *t, *mod, *invmod;

            if (fmpz_get_si(val) > 0)
            {
                mod = _fmpz_vec_init(d + 1);
                invmod = _fmpz_vec_init(d - 1);
            }
            else
            {
                mod = ctx->mod->coeffs;
                invmod = ctx->invmod->coeffs;
            }

            _fmpz_vec_scalar_fdiv_r_2exp(mod, ctx->mod->coeffs, d + 1, N - fmpz_get_si(val));
            _fmpz_vec_scalar_fdiv_r_2exp(invmod, ctx->invmod->coeffs, d - 1, N - fmpz_get_si(val));

            if (x == y)
            {
                t = _fmpz_vec_init(2 * d - 1);
            }
            else
            {
                padic_poly_fit_length(x, 2 * d - 1);
                t = x->coeffs;
            }

            _qadic_dense_pow_char_2(t, y->coeffs, y->length, e,
                                    mod, invmod, d + 1, N - fmpz_get_si(val));
            x->val = fmpz_get_si(val);

            if (x == y)
            {
                _fmpz_vec_clear(x->coeffs, x->alloc);
                x->coeffs = t;
                x->alloc  = 2 * d - 1;
                x->length = d;
            }
            else
            {
                _padic_poly_set_length(x, d);
            }
            _padic_poly_normalise(x);

            if (fmpz_get_si(val) > 0)
            {
                _fmpz_vec_clear(mod, d + 1);
                _fmpz_vec_clear(invmod, d - 1);
            }
        }
        else 
        {
            const long d = qadic_dense_ctx_degree(ctx);
            fmpz *t, *mod, *invmod;
            fmpz_t pow;
            int alloc;

            alloc = _padic_ctx_pow_ui(pow, N - fmpz_get_si(val), &ctx->pctx);

            if (fmpz_get_si(val) > 0)
            {
                mod = _fmpz_vec_init(d + 1);
                invmod = _fmpz_vec_init(d - 1);

                _fmpz_vec_scalar_mod_fmpz(mod, ctx->mod->coeffs, d + 1, pow);
                _fmpz_vec_scalar_mod_fmpz(invmod, ctx->invmod->coeffs, d - 1, pow);
            }
            else
            {
                mod = ctx->mod->coeffs;
                invmod = ctx->invmod->coeffs;
            }

            if (x == y)
            {
                t = _fmpz_vec_init(2 * d - 1);
            }
            else
            {
                padic_poly_fit_length(x, 2 * d - 1);
                t = x->coeffs;
            }

            _qadic_dense_pow(t, y->coeffs, y->length, e,
                             mod, invmod, d + 1, pow);
            x->val = fmpz_get_si(val);

            if (x == y)
            {
                _fmpz_vec_clear(x->coeffs, x->alloc);
                x->coeffs = t;
                x->alloc  = 2 * d - 1;
                x->length = d;
            }
            else
            {
                _padic_poly_set_length(x, d);
            }
            _padic_poly_normalise(x);

            if (fmpz_get_si(val) > 0)
            {
                _fmpz_vec_clear(mod, d + 1);
                _fmpz_vec_clear(invmod, d - 1);
            }

            if (alloc)
                fmpz_clear(pow);
        }
        fmpz_clear(val);
    }
}

