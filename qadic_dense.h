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

#ifndef QADIC_DENSE_H
#define QADIC_DENSE_H

#undef ulong /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long

#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"
#include "padic.h"
#include "padic_poly.h"

/* Data types and context ****************************************************/

#define qadic_dense_t padic_poly_t

#define qadic_dense_struct padic_poly_struct

#define qadic_dense_val(op) ((op)->val)

typedef struct
{
    padic_ctx_struct pctx;

    padic_poly_t mod;

    padic_poly_t invmod;

    char *var;
}
qadic_dense_ctx_struct;

typedef qadic_dense_ctx_struct qadic_dense_ctx_t[1];

int fmpz_poly_init_conway(fmpz_poly_t poly,
                          const fmpz_t p, long d);

int padic_poly_init_conway(padic_poly_t poly,
                           const fmpz_t p, long d,
                           const padic_ctx_t ctx);

void qadic_dense_ctx_init_conway(qadic_dense_ctx_t ctx, 
                           const fmpz_t p, long d, long N, const char *var, 
                           enum padic_print_mode mode);

void qadic_dense_ctx_clear(qadic_dense_ctx_t ctx);

void _qadic_dense_ctx_init_inv(fmpz *invmod, const fmpz *mod,
                               const padic_ctx_t pctx, long d, long N);

void qadic_dense_ctx_init_inv(qadic_dense_ctx_t ctx,
                              const padic_ctx_t pctx, long d, long N);

static __inline__ long qadic_dense_ctx_degree(const qadic_dense_ctx_t ctx)
{
    return padic_poly_degree(ctx->mod);
}

static __inline__ void qadic_dense_ctx_print(const qadic_dense_ctx_t ctx)
{
    printf("p    = "), fmpz_print((&ctx->pctx)->p), printf("\n");
    printf("d    = %ld\n", padic_poly_degree(ctx->mod));
    printf("N    = %ld\n", (&ctx->pctx)->N);
    printf("f(X) = ");
    padic_poly_print(ctx->mod, &ctx->pctx);
    printf("\n");
}

/* Memory management *********************************************************/

static __inline__ void qadic_dense_init(qadic_dense_t x)
{
    padic_poly_init(x);
}

static __inline__ void qadic_dense_clear(qadic_dense_t x)
{
    padic_poly_clear(x);
}

static __inline__ void
_qadic_dense_reduce_no_mod_rem(fmpz* R, const fmpz *A, long lenA,
                               const fmpz *B, const fmpz *Binv, long lenB)
{
    /* FMPZ_VEC_NORM(A, lenA); */

    if (lenA >= lenB)
    {
        _fmpz_poly_rem(R, A, lenA, B, lenB);
    }
    else
    {
        _fmpz_vec_set(R, A, lenA);
    }
}

static __inline__ void 
_qadic_dense_reduce_rem(fmpz* R, fmpz *A, long lenA,
                        const fmpz *B, const fmpz *Binv, long lenB,
                        const fmpz_t pN)
{
    /* FMPZ_VEC_NORM(A, lenA); */

    if (lenA >= lenB)
    {
        /* 1L should be replaced by a proper inverse of the leading coefficient... */
        const fmpz_t p1 = {1L};
        _fmpz_mod_poly_rem(R, A, lenA, B, lenB, p1, pN);
        /*_qadic_dense_reduce_no_mod(R, A, lenA, B, lenB);
          _fmpz_vec_scalar_mod_fmpz(R, R, lenA, pN);*/
    }
    else
    {
        _fmpz_vec_scalar_mod_fmpz(R, A, lenA, pN);
    }
}

static __inline__ void
_qadic_dense_reduce_no_mod(fmpz* R, const fmpz *A, long lenA,
                           const fmpz *B, const fmpz *Binv, long lenB)
{
    /* FMPZ_VEC_NORM(A, lenA); */

    if (lenA >= lenB)
    {
        const long m = lenA - lenB;
        long i;
        fmpz *Arev, *Q;

        Q = _fmpz_vec_init(m + 1);

        Arev = flint_malloc(lenA * sizeof(fmpz));
        for (i = 0; i < lenA; i++)
        {
            Arev[i] = A[lenA - 1 - i];
        }

        _fmpz_poly_mullow(Q, Arev, lenA, Binv, lenB, m + 1);
        _fmpz_poly_reverse(Q, Q, m + 1, m + 1);

        if (m < lenB)
            _fmpz_poly_mullow(R, B, lenB, Q, m + 1, lenB);
        else
            _fmpz_poly_mullow(R, Q, m + 1, B, lenB, lenB);
        _fmpz_poly_sub(R, A, lenA, R, lenB);

        _fmpz_vec_clear(Q, m + 1);

        flint_free(Arev);
    }
    else
    {
        _fmpz_vec_set(R, A, lenA);
    }
}

static __inline__ void 
_qadic_dense_reduce(fmpz* R, fmpz *A, long lenA,
                    const fmpz *B, const fmpz *Binv, long lenB,
                    const fmpz_t pN)
{
    /* FMPZ_VEC_NORM(A, lenA); */

    if (lenA >= lenB)
    {
        const long m = lenA - lenB;
        long i;
        fmpz *Arev, *Q;

        Q = _fmpz_vec_init(m + 1);
        Arev = flint_malloc(lenA * sizeof(fmpz));
        for (i = 0; i < lenA; i++)
        {
            Arev[i] = A[lenA - 1 - i];
        }

        _fmpz_poly_mullow(Q, Arev, lenA, Binv, lenB - 2, m + 1);
        _fmpz_vec_scalar_mod_fmpz(Q, Q, m + 1, pN);
        _fmpz_poly_reverse(Q, Q, m + 1, m + 1);

        if (lenB >= m + 1)
            _fmpz_poly_mullow(R, B, lenB, Q, m + 1, lenB);
        else
            _fmpz_poly_mullow(R, Q, m + 1, B, lenB, lenB);
        _fmpz_poly_sub(R, A, lenA, R, lenB);

        _fmpz_vec_scalar_mod_fmpz(R, R, lenB, pN);
        _fmpz_vec_zero(R + lenB, m);

        _fmpz_vec_clear(Q, m + 1);

        flint_free(Arev);
    }
    else
    {
        _fmpz_vec_scalar_mod_fmpz(R, A, lenA, pN);
    }
}

static __inline__ void 
_qadic_dense_reduce_char_2(fmpz* R, fmpz *A, long lenA,
                    const fmpz *B, const fmpz *Binv, long lenB,
                    long N)
{
    /* FMPZ_VEC_NORM(A, lenA); */

    if (lenA >= lenB)
    {
        const long m = lenA - lenB;
        long i;
        fmpz *Arev, *Q;

        Q = _fmpz_vec_init(m + 1);

        Arev = flint_malloc(lenA * sizeof(fmpz));
        for (i = 0; i < lenA; i++)
        {
            Arev[i] = A[lenA - 1 - i];
        }

        _fmpz_poly_mullow(Q, Arev, lenA, Binv, lenB - 2, m + 1);
        _fmpz_vec_scalar_fdiv_r_2exp(Q, Q, m + 1, N);
        _fmpz_poly_reverse(Q, Q, m + 1, m + 1);

        if (lenB >= m + 1)
            _fmpz_poly_mullow(R, B, lenB, Q, m + 1, lenB);
        else
            _fmpz_poly_mullow(R, Q, m + 1, B, lenB, lenB);
        _fmpz_poly_sub(R, A, lenA, R, lenB);

        _fmpz_vec_scalar_fdiv_r_2exp(R, R, lenB, N);

        _fmpz_vec_clear(Q, m + 1);

        flint_free(Arev);
    }
    else
    {
        _fmpz_vec_scalar_fdiv_r_2exp(R, A, lenA, N);
    }
}

static __inline__ void qadic_dense_reduce(qadic_dense_t x, const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;
    const long d = qadic_dense_ctx_degree(ctx);

    if (x->length == 0 || x->val >= N)
    {
        padic_poly_zero(x);
    }
    else if (!COEFF_IS_MPZ(*(&ctx->pctx)->p) && *(&ctx->pctx)->p == 2L)
    {
        fmpz *t, *mod, *invmod;

        t      = _fmpz_vec_init(x->length);

        if (x->val > 0)
        {
            mod    = _fmpz_vec_init(d + 1);
            invmod = _fmpz_vec_init(d - 1);

            _fmpz_vec_scalar_fdiv_r_2exp(mod, ctx->mod->coeffs, d + 1, N - x->val);
            _fmpz_vec_scalar_fdiv_r_2exp(invmod, ctx->invmod->coeffs, d - 1, N - x->val);
        }
        else
        {
            mod = ctx->mod->coeffs;
            invmod = ctx->invmod->coeffs;
        }

        _qadic_dense_reduce_char_2(t, x->coeffs, x->length,
                            mod, invmod, d + 1,
                            N - x->val);
        _fmpz_vec_set(x->coeffs, t, FLINT_MIN(x->length, d));
        _padic_poly_set_length(x, FLINT_MIN(x->length, d));
        _padic_poly_normalise(x);
        padic_poly_canonicalise(x, (&ctx->pctx)->p);

        _fmpz_vec_clear(t, x->length);
        if (x->val > 0)
        {
            _fmpz_vec_clear(mod, d + 1);
            _fmpz_vec_clear(invmod, d - 1);
        }
    }
    else
    {
        fmpz_t pow;
        fmpz *t, *mod, *invmod;
        int alloc;

        t      = _fmpz_vec_init(x->length);

        alloc = _padic_ctx_pow_ui(pow, (&ctx->pctx)->N - x->val, &ctx->pctx);

        if (x->val)
        {
            mod    = _fmpz_vec_init(d + 1);
            invmod = _fmpz_vec_init(d - 1);

            _fmpz_vec_scalar_mod_fmpz(mod, ctx->mod->coeffs, d + 1, pow);
            _fmpz_vec_scalar_mod_fmpz(invmod, ctx->invmod->coeffs, d - 1, pow);
        }
        else
        {
            mod = ctx->mod->coeffs;
            invmod = ctx->invmod->coeffs;
        }

        _qadic_dense_reduce(t, x->coeffs, x->length,
                            mod, invmod, d + 1,
                            pow);
        _fmpz_vec_set(x->coeffs, t, FLINT_MIN(x->length, d));
        _padic_poly_set_length(x, FLINT_MIN(x->length, d));
        _padic_poly_normalise(x);
        padic_poly_canonicalise(x, (&ctx->pctx)->p);

        if (alloc)
            fmpz_clear(pow);

        _fmpz_vec_clear(t, x->length);
        if (x->val > 0)
        {
            _fmpz_vec_clear(mod, d + 1);
            _fmpz_vec_clear(invmod, d - 1);
        }
    }
}

void qadic_dense_scalar_mod_ppow(qadic_dense_t rop, const qadic_dense_t op, long e, 
                           const qadic_dense_ctx_t ctx);

/* Randomisation *************************************************************/

static __inline__ void 
qadic_dense_randtest(qadic_dense_t x, flint_rand_t state, const qadic_dense_ctx_t ctx)
{
    padic_poly_randtest(x, state, qadic_dense_ctx_degree(ctx), &ctx->pctx);
}

static __inline__ void 
qadic_dense_randtest_not_zero(qadic_dense_t x, flint_rand_t state, const qadic_dense_ctx_t ctx)
{
    padic_poly_randtest_not_zero(x, state, qadic_dense_ctx_degree(ctx), 
                                 &ctx->pctx);
}

static __inline__ void 
qadic_dense_randtest_val(qadic_dense_t x, flint_rand_t state, long val, 
                   const qadic_dense_ctx_t ctx)
{
    padic_poly_randtest_val(x, state, val, qadic_dense_ctx_degree(ctx), &ctx->pctx);
}

static __inline__ void 
qadic_dense_randtest_int(qadic_dense_t x, flint_rand_t state, const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (N <= 0)
    {
        padic_poly_zero(x);
    }
    else
    {
        padic_poly_randtest_val(x, state, n_randint(state, N), 
                                qadic_dense_ctx_degree(ctx), &ctx->pctx);
    }
}

/* Assignments and conversions ***********************************************/

static __inline__ void qadic_dense_zero(qadic_dense_t x)
{
    padic_poly_zero(x);
}

static __inline__ void qadic_dense_one(qadic_dense_t x, const qadic_dense_ctx_t ctx)
{
    padic_poly_one(x, &ctx->pctx);
}

static __inline__ void qadic_dense_gen(qadic_dense_t x, const qadic_dense_ctx_t ctx)
{
    const long d = qadic_dense_ctx_degree(ctx);

    if (d > 1)
    {
        if ((&ctx->pctx)->N > 0)
        {
            padic_poly_fit_length(x, 2);
            fmpz_zero(x->coeffs + 0);
            fmpz_one(x->coeffs + 1);
            _padic_poly_set_length(x, 2);
            x->val = 0;
        }
        else
        {
            padic_poly_zero(x);
        }
    }
    else
    {
        printf("Exception (qadic_dense_gen).  Extension degree d = 1.\n");
        abort();
    }
}

static __inline__ 
void qadic_dense_set_ui(qadic_dense_t rop, ulong op, const qadic_dense_ctx_t ctx)
{
    padic_poly_set_ui(rop, op, &ctx->pctx);
}

static __inline__ int 
qadic_dense_get_padic(padic_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx)
{
    if (op->length > 0)
    {
        if (_fmpz_vec_is_zero(op->coeffs + 1, op->length - 1))
        {
            fmpz_set(padic_unit(rop), op->coeffs + 0);
            padic_val(rop) = op->val;
            _padic_canonicalise(rop, &ctx->pctx);
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        padic_zero(rop);
        return 1;
    }
}

static __inline__ void qadic_dense_set(qadic_dense_t x, const qadic_dense_t y)
{
    padic_poly_set(x, y);
}

void qadic_dense_set_fmpz_poly(qadic_dense_t rop, const fmpz_poly_t op, 
                         const qadic_dense_ctx_t ctx);

/* Comparison ****************************************************************/

static __inline__ int 
qadic_dense_is_zero(const qadic_dense_t x)
{
    return padic_poly_is_zero(x);
}

static __inline__ int 
qadic_dense_is_one(const qadic_dense_t x, const qadic_dense_ctx_t ctx)
{
    return padic_poly_is_one(x, &ctx->pctx);
}

static __inline__ int 
qadic_dense_equal(const qadic_dense_t x, const qadic_dense_t y)
{
    return padic_poly_equal(x, y);
}

/* Basic arithmetic **********************************************************/

static __inline__ void 
qadic_dense_add(qadic_dense_t x, const qadic_dense_t y, const qadic_dense_t z, const qadic_dense_ctx_t ctx)
{
    padic_poly_add(x, y, z, &ctx->pctx);
}

static __inline__ void 
qadic_dense_sub(qadic_dense_t x, const qadic_dense_t y, const qadic_dense_t z, const qadic_dense_ctx_t ctx)
{
    padic_poly_sub(x, y, z, &ctx->pctx);
}

static __inline__ void 
qadic_dense_neg(qadic_dense_t x, const qadic_dense_t y, const qadic_dense_ctx_t ctx)
{
    padic_poly_neg(x, y, &ctx->pctx);
}

void qadic_dense_mul(qadic_dense_t x, const qadic_dense_t y, const qadic_dense_t z,
                     const qadic_dense_ctx_t ctx);

void _qadic_dense_inv(fmpz *rop, const fmpz *op, long len, 
                const fmpz *mod, const fmpz *invmod, long lenmod, 
                const fmpz_t p, long N);

void qadic_dense_inv(qadic_dense_t x, const qadic_dense_t y, const qadic_dense_ctx_t ctx);

void qadic_dense_sqr(qadic_dense_t x, const qadic_dense_t y, const qadic_dense_ctx_t ctx);

void _qadic_dense_pow(fmpz *rop, const fmpz *op, long len, const fmpz_t e, 
                const fmpz *mod, const fmpz *invmod, long lenmod, 
                const fmpz_t p);

void qadic_dense_pow(qadic_dense_t x, const qadic_dense_t y, const fmpz_t e, const qadic_dense_ctx_t ctx);

/* Special functions *********************************************************/

void _qadic_dense_exp_rectangular(fmpz *rop, const fmpz *op, long v, long len, 
                            const fmpz *mod, const fmpz *invmod, long lenmod, 
                            const fmpz_t p, long N, const fmpz_t pN);

int qadic_dense_exp_rectangular(qadic_dense_t rop, const qadic_dense_t op, 
                          const qadic_dense_ctx_t ctx);

void _qadic_dense_exp_balanced(fmpz *rop, const fmpz *op, long v, long len, 
                            const fmpz *mod, const fmpz *invmod, long lenmod, 
                            const fmpz_t p, long N, const fmpz_t pN);

int qadic_dense_exp_balanced(qadic_dense_t rop, const qadic_dense_t op, 
                          const qadic_dense_ctx_t ctx);

void _qadic_dense_exp(fmpz *rop, const fmpz *op, long v, long len, 
                const fmpz *mod, const fmpz *invmod, long lenmod, 
                const fmpz_t p, long N, const fmpz_t pN);

int qadic_dense_exp(qadic_dense_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx);

void _qadic_dense_log_rectangular(fmpz *z, const fmpz *y, long v, long len, 
                            const fmpz *mod, const fmpz *invmod, long lenmod, 
                            const fmpz_t p, long N, const fmpz_t pN);

int qadic_dense_log_rectangular(qadic_dense_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx);

void _qadic_dense_log_balanced(fmpz *z, const fmpz *y, long len, 
                         const fmpz *mod, const fmpz *invmod, long lenmod, 
                         const fmpz_t p, long N, const fmpz_t pN);

int qadic_dense_log_balanced(qadic_dense_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx);

void _qadic_dense_log(fmpz *z, const fmpz* y, long v, long len,
                      const fmpz* mod, const fmpz *invmod, long lenmod,
                      const fmpz_t p, long N, const fmpz_t pN);

int qadic_dense_log(qadic_dense_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx);

void qadic_dense_frobenius(qadic_dense_t rop, const qadic_dense_t op, long e, const qadic_dense_ctx_t ctx);

void _qadic_dense_teichmuller(fmpz *rop, const fmpz *op, long len, 
                        const fmpz *mod, const fmpz *invmod, long lenmod, 
                        const fmpz_t p, long N);

void qadic_dense_teichmuller(qadic_dense_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx);

void _qadic_dense_trace(fmpz_t rop, const fmpz *op, long len, 
                  const fmpz *mod, const fmpz *invmod, long lenmod, 
                  const fmpz_t pN);

void qadic_dense_trace(padic_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx);

void _qadic_dense_norm(fmpz_t rop, const fmpz *op, long len, 
                 const fmpz *mod, const fmpz *invmod, long lenmod, 
                 const fmpz_t p, long N);

void qadic_dense_norm(padic_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx);

void _qadic_dense_norm_analytic(fmpz_t rop, const fmpz *y, long v, long len, 
                          const fmpz *mod, const fmpz *invmod, long lenmod, 
                          const fmpz_t p, long N);

void qadic_dense_norm_analytic(padic_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx);

void _qadic_dense_norm_resultant(fmpz_t rop, const fmpz *op, long len, 
                           const fmpz *mod, const fmpz *invmod, long lenmod, 
                           const fmpz_t p, long N);

void qadic_dense_norm_resultant(padic_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx);

int qadic_dense_sqrt(qadic_dense_t rop, const qadic_dense_t op, const qadic_dense_ctx_t ctx);

/* Output ********************************************************************/

int qadic_dense_fprint_pretty(FILE *file, const qadic_dense_t op, const qadic_dense_ctx_t ctx);

static __inline__ int 
qadic_dense_print_pretty(const qadic_dense_t op, const qadic_dense_ctx_t ctx)
{
    return qadic_dense_fprint_pretty(stdout, op, ctx);
}

#endif

