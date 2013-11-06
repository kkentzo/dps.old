#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "sampling.h"
#include "clib.h"


//===============================================================================
//                       UTILITY FUNCTIONS
//===============================================================================



// mutate the supplied value by sampling from a uniform distribution in [VAL-DELTA, VAL+DELTA]
// make sure that the returned value is within the supplied LO and HI and boundaries
// PDC indicates whether there exist periodic boundary conditions
double mutate(double val, double delta, double lo, double hi, int pbc, void *rng) {

    double x1, x2, new_val;

    if (pbc) {

	x1 = val - delta;
	x2 = val + delta;
	
    } else {

	x1 = MYMAX(val-delta, lo);
	x2 = MYMIN(val+delta, hi);

    }

    new_val = x1 + (x2 - x1) * gsl_rng_uniform(rng);

    if (pbc) {
	if (new_val < lo)
	    new_val += (hi-lo);
	else if (new_val > hi)
	    new_val -= (hi-lo);
    }

    return new_val;
    
}





// normalizes the vector V (of length N) so that \sum_i v_i = 1
// if VSUM > 0, it is used as the vector's sum 
// otherwise the sum of the vector is computed
void vnormalize(double *v, int n, double vsum) {

    int i;

    // compute sum??
    if (vsum <= 0) 
	for (i=0, vsum=0; i<n; vsum += v[i++]) ;

    // normalize vector
    for (i=0; i<n; v[i++] /= vsum) ;
    
}






// shuffles (randomizes the order of) the contents of the VEC
// algorithm from :http://benpfaff.org/writings/clc/shuffle.html

// see also http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
void shuffle(void *v, int n, int size, void *rng) {

    if (n <= 1)
	return;

    int i, j;
    unsigned long rng_max = gsl_rng_max((gsl_rng *)rng);

    void *buf = malloc(size);

    for (i=0; i<n-1; i++) {

	// INVESTIGATE WHETHER SOMETHING LIKE THIS WORKS AS WELL
	//j = i + (n-i) * gsl_rng_uniform_int(rng, n);

	j = i + gsl_rng_get(rng) / (rng_max / (n - i) + 1);

	if (i==j)
	    continue;

	// swap i^th and j^th items
	memcpy(buf, v+i*size, size);
	memcpy(v+i*size, v+j*size, size);
	memcpy(v+j*size, buf, size);
	
    }

    free(buf);
    
}






//=================================
//     SAMPLING - NO EXCLUSION
//=================================



// equal sampling with replacement
static void sample_replace(int k, int n, int *dest, void *rng) {

    int i;
    for (i=0; i<k; i++)
	dest[i] = n * gsl_rng_uniform(rng);
    
}




// equal sampling without replacement (i.e. with removal)
static void sample_noreplace(int k, int n, int *dest, void *rng) {

    int i, j;
    int *x = calloc(n, sizeof(int));

    for (i = 0; i < n; i++)
	x[i] = i;
    
    for (i = 0; i < k; i++) {
	j = n * gsl_rng_uniform(rng);
	dest[i] = x[j];
	x[j] = x[--n];
    }

    free(x);
    
}





// unequal sampling with replacement
static void sample_unequal_replace(int k, int n, int *dest, double *src, void *rng) {


    /* Create the alias tables.
       The idea is that for HL[0] ... L-1 label the entries with q < 1
       and L ... H[n-1] label those >= 1.
       By rounding error we could have q[i] < 1. or > 1. for all entries.
    */

    int *a = calloc(n, sizeof(int));
    int *HL = calloc(n, sizeof(int));
    double *q = calloc(n, sizeof(double));

    int *H, *L;
    double rU;
    int i, j, m;
    

    H = HL - 1;
    L = HL + n;
    for (i = 0; i < n; i++) {
	q[i] = src[i] * n;
	if (q[i] < 1.)
	    *++H = i;
	else
	    *--L = i;
    }
    if (H >= HL && L < HL + n) { /* So some q[i] are >= 1 and some < 1 */
	for (m = 0; m < n - 1; m++) {
	    i = HL[m];
	    j = *L;
	    a[i] = j;
	    q[j] += q[i] - 1;
	    if (q[j] < 1.) L++;
	    if(L >= HL + n) break; /* now all are >= 1 */
	}
    }
    
    for (i = 0; i < n; i++)
	q[i] += i;

    /* generate sample */
    for (i = 0; i < k; i++) {
	rU = gsl_rng_uniform(rng) * n;
	m = (int) rU;
	dest[i] = (rU < q[m]) ? m : a[m];
    }

    // free arrays
    free(a);
    free(HL);
    free(q);
    
}




// unequal sampling without replacement (i.e. with removal)
static void sample_unequal_noreplace(int k, int n, int *dest, double *src, void *rng) {

    double rT, mass, totalmass;
    int i, j, m, n1;

    int *ind = calloc(n, sizeof(int));

    /* Record element identities */
    for (i = 0; i < n; i++)
    	ind[i] = i;

    /* Sort probabilities into descending order */
    /* Order element identities in parallel */
    qsort_w_indices(src, ind, n, sizeof(double), dcmp_desc);

    /* Compute the sample */
    totalmass = 1;
    for (i = 0, n1 = n-1; i < k; i++, n1--) {
    	rT = totalmass * gsl_rng_uniform(rng);
    	mass = 0;
    	for (j = 0; j < n1; j++) {
    	    mass += src[j];
    	    if (rT <= mass)
    		break;
    	}
    	dest[i] = ind[j];
    	totalmass -= src[j];
    	for(m = j; m < n1; m++) {
    	    src[m] = src[m + 1];
    	    ind[m] = ind[m + 1];
    	}
    }

    free(ind);
    
    
}






// sample k integers from [0,n) and store them in dest
// if src==NULL, all options are considered with equal probability
// otherwise, each option i is considered with probability src[i]
// As such it must be that len(dest) == k and len(src) == n
// REPLACE determines whether sampling occurs with replacement or not
// NORM determines whether src needs to be normalized (so that \sum src_i = 1)
void sample(int k, int n, int *dest, double *src, int replace, int norm, void *rng) {

    int i;

    // we can't sample more than n without replacement
    assert(! (!replace && k > n));

    // take the shortcut when k==n with no replacement
    // just copy all the indices over
    if (!replace && k == n) {
	// select and return all numbers
	for (i=0; i<n; i++)
	    dest[i] = i;
	return;
    }
    
    // normalize probability vector??
    if (src && norm) 
	vnormalize(src, n, -1);

    if (src) {
	if (replace)
	    // unequal sampling with replacement
	    sample_unequal_replace(k, n, dest, src, rng);
	else
	    // unequal sampling without replacement
	    sample_unequal_noreplace(k, n, dest, src, rng);
    } else {

	if (replace)
	    // equal sampling with replacement
	    sample_replace(k, n, dest, rng);
	else
	    // equal sampling without replacement
	    sample_noreplace(k, n, dest, rng);
	
    }
    
}




//=================================
//   SAMPLING - WITH EXCLUSION
//=================================

// sample k integers from [0,n) and store them in dest
// excluding the values stored in ev (with length en)
// all options are considered with equal probability
void sample_excluding(int k, int n, int *dest, int *ev, int en, void *rng) {

    int rnd;
    int i, j;

    for (i=0; i<k; i++) {

	rnd = gsl_rng_uniform_int((gsl_rng *)rng, n - en);
	
	for (j=0; j<en; j++) {
	    if (rnd < ev[i])
		break;
	    rnd ++;
	}

	dest[i] = rnd;

    }
    
}



