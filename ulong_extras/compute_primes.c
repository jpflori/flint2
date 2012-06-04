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

    Copyright (C) 2009 Tom Boothby
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#undef ulong /* prevent clash with standard library */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>
#define ulong unsigned long
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"

unsigned int * sieve;

const unsigned int flint_primes_small[] =
{
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
    101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
    191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,
    281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,
    389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,
    491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,
    607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,
    719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,
    829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,
    953,967,971,977,983,991,997,1009,1013,1019,1021
};

mp_limb_t * flint_primes;
mp_limb_t flint_primes_cutoff = 0;

double * flint_prime_inverses;

ulong flint_num_primes = 0;
pthread_mutex_t flint_num_primes_mutex;

void n_compute_primes(ulong num)
{
    ulong num_primes, primes_cutoff;
    ulong sieve_size;
    unsigned int i, j, p, q, oldq = 0;

    if (flint_num_primes >= num) return;

    pthread_mutex_lock(&flint_num_primes_mutex);
    if (flint_num_primes >= num) /* someone may have changed this before we locked */
    {
        pthread_mutex_unlock(&flint_num_primes_mutex);
        return; 
    }

    num_primes = FLINT_MAX(flint_num_primes, FLINT_NUM_PRIMES_SMALL);

    if (!flint_num_primes) flint_primes_cutoff = FLINT_PRIMES_SMALL_CUTOFF;

    num = FLINT_MAX(num, 16384);
    n_nth_prime_bounds(&primes_cutoff, &primes_cutoff, num+1);
    /* bad things will happen if primes_cutoff equals a prime */
    primes_cutoff += primes_cutoff & 1;

    if (primes_cutoff > FLINT_PRIMES_SMALL_CUTOFF*FLINT_PRIMES_SMALL_CUTOFF)
    {
        printf("Exception: cannot precompute sufficiently many primes!\n");
        abort();
    }

    sieve_size = primes_cutoff/2 - flint_primes_cutoff/2;
    sieve = (unsigned int *) flint_malloc(sizeof(unsigned int)*sieve_size);
   
    for (j = 0; j < sieve_size; j++)
      sieve[j] = 1;

    for (i = 1; ; i++) /* skip 2 */
    {
        p = flint_primes_small[i];
        if (p*p > flint_primes_cutoff) break;

        q = flint_primes_cutoff - (flint_primes_cutoff%p); /* first multiple of p before sieve start */
        q = q + (1 + q%2)*p; /* first odd multiple of p after sieve start */
        q = q/2 - flint_primes_cutoff/2;   /* index of first odd multiple of p left in the sieve */

        while (q < sieve_size)
        {
            sieve[q] = 0;
            q += p;
        }
    }

    for ( ; (i < FLINT_NUM_PRIMES_SMALL); i++)
    {
        p = flint_primes_small[i];
        if (p*p > primes_cutoff) break;

        q = (p*p)/2 - flint_primes_cutoff/2;   /* index of first multiple of p left in the sieve */

        for (j = oldq; j < q; j++) /* count new primes up to p^2 */
        {
            if (sieve[j]) num_primes++;
        }

        oldq = q;

        while (q < sieve_size)
        {
            sieve[q] = 0;
            q += p;
        }
    }

    for (j = oldq; j < sieve_size; j++) /* there may be some primes after the last p^2 */
    {
        if (sieve[j]) num_primes++;
    }

    if (!flint_num_primes) 
    {
        flint_primes = (mp_limb_t *) flint_malloc(sizeof(mp_limb_t)*num_primes);
        flint_prime_inverses = (double *) flint_malloc(sizeof(double)*num_primes);
    }
    else
    {
        flint_primes = flint_realloc(flint_primes, sizeof(mp_limb_t)*num_primes);
        flint_prime_inverses = flint_realloc(flint_prime_inverses, sizeof(double)*num_primes);
    }

    if (!flint_num_primes)
    {
        for (i = 0; i < FLINT_NUM_PRIMES_SMALL; i++) /* write out small primes */
        {
            p = flint_primes[i] = flint_primes_small[i];
            flint_prime_inverses[i] = n_precompute_inverse(p);
            flint_num_primes = FLINT_NUM_PRIMES_SMALL;
        }
    }
   
    i = flint_num_primes;

    for (j = 0; j < sieve_size; j++) /* write out primes computed in sieve */
    {
        if (sieve[j]) 
        {
            p = flint_primes[i] = 2*(j + flint_primes_cutoff/2) + 1;
            flint_prime_inverses[i++] = n_precompute_inverse(p);
        }
    }

    flint_primes_cutoff = primes_cutoff;
    if (flint_primes_cutoff & 1) flint_primes_cutoff++;

    flint_num_primes = num_primes;

    flint_free(sieve);

    pthread_mutex_unlock(&flint_num_primes_mutex);
}
