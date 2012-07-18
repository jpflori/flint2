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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

#define E(j,k) nmod_poly_mat_entry(B,j,k)

static __inline__ void
nmod_poly_mat_swap_rows(nmod_poly_mat_t mat, long * perm, long r, long s)
{
    if (r != s)
    {
        nmod_poly_struct * u;
        long t;

        if (perm != NULL)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u; 
    }
}

long
nmod_poly_mat_rref(nmod_poly_mat_t B, nmod_poly_t den,
    const nmod_poly_mat_t A)
{
    nmod_poly_t t;
    long m, n, j, k, rank, r, pivot_row, pivot_col;
    int sign;

    if (nmod_poly_mat_is_empty(A))
    {
        nmod_poly_one(den);
        return 0;
    }

    nmod_poly_mat_set(B, A);
    m = B->r;
    n = B->c;
    rank = pivot_row = pivot_col = 0;
    sign = 1;

    nmod_poly_init(t, nmod_poly_mat_modulus(A));

    while (pivot_row < m && pivot_col < n)
    {
        r = nmod_poly_mat_find_pivot_partial(B, pivot_row, m, pivot_col);

        if (r == -1)
        {
            pivot_col++;
            continue;
        }
        else if (r != pivot_row)
        {
            nmod_poly_mat_swap_rows(B, NULL, pivot_row, r);
            sign = -sign;
        }
        rank++;

        for (j = 0; j < m; j++)
        {
            if (j == pivot_row)
                continue;

            for (k = j > pivot_row ? pivot_col + 1 : j; k < n; k++)
            {
                if (k == pivot_col)
                    continue;

                nmod_poly_mul(E(j, k), E(j, k), E(pivot_row, pivot_col));
                nmod_poly_mul(t, E(j, pivot_col), E(pivot_row, k));
                nmod_poly_sub(E(j, k), E(j, k), t);

                if (pivot_row > 0)
                    nmod_poly_div(E(j, k), E(j, k), den);
            }

            nmod_poly_zero(E(j, pivot_col));
        }

        nmod_poly_set(den, E(pivot_row, pivot_col));
        pivot_row++;
        pivot_col++;
    }

    nmod_poly_clear(t);

    if (sign < 0)
    {
        nmod_poly_neg(den, den);
        nmod_poly_mat_neg(B, B);
    }

    return rank;
}
