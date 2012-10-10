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

    Copyright (C) 2012 Jean-Pierre Flori

******************************************************************************/

#include "padic_poly.h"
#include "qadic_dense.h"

void _qadic_dense_ctx_init_inv(fmpz *invmod, const fmpz *mod,
                               const padic_ctx_t pctx, long d, long N)
{
    /* We assume for now the leading coefficient of the modulus is one */
    const fmpz_t one = {1L};
    int alloc;
    long i;
    fmpz_t pow;
    fmpz *revmod;

    alloc = _padic_ctx_pow_ui(pow, N, pctx);

    revmod = flint_malloc((d - 1) * sizeof(fmpz));
    for (i = 0; i < d - 1; i++)
    {
        revmod[i] = mod[d - i];
    }

    _fmpz_mod_poly_inv_series_newton(invmod, revmod, d - 1, one, pow);

    flint_free(revmod);

    if (alloc)
        fmpz_clear(pow);
}

void qadic_dense_ctx_init_inv(qadic_dense_ctx_t ctx,
                              const padic_ctx_t pctx, long d, long N)
{
    padic_poly_init2(ctx->invmod, d - 1);
    _qadic_dense_ctx_init_inv(ctx->invmod->coeffs, ctx->mod->coeffs, pctx, d, N);
    _padic_poly_set_length(ctx->invmod, d - 1);
    _padic_poly_normalise(ctx->invmod);
}
