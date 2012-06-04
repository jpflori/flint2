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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t 
n_powmod2_preinv(mp_limb_t a, mp_limb_signed_t exp, mp_limb_t n, mp_limb_t ninv)
{
    mp_limb_t x, y;
    mp_limb_t e;
    
    if (n == 1UL) return 0UL;

    e = (exp < 0L ? -exp : exp);

    x = 1UL;
    y = a;
    while (e) 
    {
        if (e & 1) x = n_mulmod2_preinv(x, y, n, ninv);
        e = e >> 1;
        if (e) y = n_mulmod2_preinv(y, y, n, ninv);
    }

    return (exp < 0L ? n_invmod(x, n) : x);
} 
