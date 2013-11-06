#include <string.h>
#include <math.h>
#include <assert.h>

#include "clib.h"

#include "nvar.h"




nvar_t *nvar_init(dim) {

    nvar_t *nvar = malloc(sizeof(nvar_t));

    nvar->dim = dim;

    nvar->n = 0;

    nvar->wsums = calloc(dim, sizeof(double));

    nvar->sum_wprods = sm_new(dim, 0);

    // allocate space for the array of pointers to variable labels
    nvar->labels = calloc(nvar->dim, sizeof(char *));

    // allocate space for the array of pointers to covariance labels
    nvar->plabels = calloc(nvar_count_covariances(nvar), sizeof(char *));

    return nvar;
    
}

void nvar_free(nvar_t *nvar) {

    // free individual labels here
    int i;
    for (i=0; i<nvar->dim; i++)
	if (nvar->labels[i])
	    free(nvar->labels[i]);

    int nlabels = nvar_count_covariances(nvar);
    for (i=0; i<nlabels; i++)
	if (nvar->plabels[i])
	    free(nvar->plabels[i]);

    
    free(nvar->wsums);

    sm_free(nvar->sum_wprods);

    // free the array 
    free(nvar->labels);
    free(nvar->plabels);

    free(nvar);

}



void nvar_reset(nvar_t *nvar) {

    // reset vector of weighted sums
    BFILL(nvar->wsums, nvar->dim, 0);
    //memset(nvar->wsums, 0, nvar->dim * sizeof(double));
    // reset product matrix
    sm_fill(nvar->sum_wprods, 0);
    // reset weight count
    nvar->n = 0;
    
}





// add variable labels (copy the strings)
void nvar_assign_labels(nvar_t *nvar, const char **labels) {

    int i, j, cnt;
    int lengths[nvar->dim];

    
    
    // figure out string lengths
    for (i=0; i<nvar->dim; i++)
	lengths[i] = strlen(labels[i]);

    // create variable labels
    for (i=0, cnt=0; i<nvar->dim; i++) {
	// reserve space for the string AND the \0 char
	nvar->labels[i] = malloc((lengths[i]+1) * sizeof(char));
	// copy string buffer
	memmove(nvar->labels[i], labels[i], lengths[i]+1);

	// create covariance labels
	for (j=i+1; j<nvar->dim; j++, cnt++) {
	    nvar->plabels[cnt] = malloc(lengths[i]+lengths[j]+2 * sizeof(char));
	    // copy string1 buffer
	    memmove(nvar->plabels[cnt], labels[i], lengths[i]);
	    // insert dot
	    nvar->plabels[cnt][lengths[i]] = '.';
	    // copy string2 buffer
	    memmove(nvar->plabels[cnt]+lengths[i]+1, labels[j], lengths[j]+1);
	}
    }
    
}




// print the labels
void nvar_print_cov_labels(nvar_t *nvar) {

    int i, j, cnt;

    for (i=0, cnt=0; i<nvar->dim; i++) 
	for (j=i+1; j<nvar->dim; j++, cnt++)
	    printf("%d %s\n", cnt, nvar->plabels[cnt]);
    
}



// add the data point specified in VALUES (of length nvar->dim) 
// with the supplied WEIGHT
void nvar_add(nvar_t *nvar, double *values, double weight) {

    int i, j;

    for (i=0; i<nvar->dim; i++) {

	// update weighted sum
	nvar->wsums[i] += weight * values[i];
	// update products
	for (j=i; j<nvar->dim; j++)
	    ((sm_t *)nvar->sum_wprods)->data[i][j-i] += weight * values[i] * values[j];
	    /* sm_set(nvar->sum_wprods, i, j, */
	    /* 	   sm_get(nvar->sum_wprods, i, j) + weight * values[i] * values[j]); */
    }

    // update weight counter
    nvar->n += weight;
    
}




// add the supplied stats (NVAR2) to NVAR1
void nvar_add_stats(nvar_t *nvar1, nvar_t *nvar2) {

    // make sure we've got the same dimensionality
    assert(nvar1->dim == nvar2->dim);

    int i, j;

    for (i=0; i<nvar1->dim; i++) {

	// update weighted sum
	nvar1->wsums[i] += nvar2->wsums[i];
	// update products
	for (j=i; j<nvar1->dim; j++)
	    ((sm_t *)nvar1->sum_wprods)->data[i][j-i] +=
		((sm_t *)nvar2->sum_wprods)->data[i][j-i];
	    /* sm_set(nvar1->sum_wprods, i, j, */
	    /* 	   sm_get(nvar1->sum_wprods, i, j) + */
	    /* 	   sm_get(nvar2->sum_wprods, i, j)); */
	
    }

    // update the numbr of observations (weight counter)
    nvar1->n += nvar2->n;
    
}




// remove the data point specified in VALUES (of length nvar->dim) 
// with the supplied WEIGHT
void nvar_remove(nvar_t *nvar, double *values, double weight) {

    int i, j;

    for (i=0; i<nvar->dim; i++) {

	// update weighted sum
	nvar->wsums[i] -= weight * values[i];
	// update products
	for (j=i; j<nvar->dim; j++)
	    ((sm_t *)nvar->sum_wprods)->data[i][j-i] -= weight * values[i] * values[j];
	    /* sm_set(nvar->sum_wprods, i, j, */
	    /* 	   sm_get(nvar->sum_wprods, i, j) - weight * values[i] * values[j]); */
    }

    // update weight counter
    nvar->n -= weight;
    
}




// calculate the mean of the specified variable
double nvar_calc_mean(nvar_t *nvar, int idx) {

    return nvar->wsums[idx] / nvar->n;
    
}





// calculate the covariance between variables IDX1 and IDX2
// if IDX1 == IDX2 returns the variable's variance 
// if (NORM) then returns the correlation coefficient between IDX1, IDX2
//            and if IDX1==IDX2 then the result is 1
double nvar_calc_covariance(nvar_t *nvar, int idx1, int idx2, int norm) {

    if (idx1 > idx2) {
	int tmp = idx1; idx1=idx2; idx2=tmp;
    }
	

    double E_x_y = ((sm_t *)nvar->sum_wprods)->data[idx1][idx2-idx1] / nvar->n;
	//sm_get(nvar->sum_wprods, idx1, idx2) / nvar->n;
    double E_x = nvar_calc_mean(nvar, idx1);
    double E_y = nvar_calc_mean(nvar, idx2);

    // === biased estimator for covariance ===
    return (E_x_y - E_x * E_y) /
    	(norm ? sqrt(nvar_calc_covariance(nvar, idx1, idx1, FALSE) *
    		     nvar_calc_covariance(nvar, idx2, idx2, FALSE)) : 1);

    // === unbiased estimator for covariance ===
    /* return nvar->n * (E_x_y - E_x * E_y) / (nvar->n - 1) / */
    /* 	(norm ? sqrt(nvar_calc_covariance(nvar, idx1, idx1, FALSE) * */
    /* 		     nvar_calc_covariance(nvar, idx2, idx2, FALSE)) : 1); */
    
}







// ================= CALCULATE MULTIPLE ELEMENT STATISTICS ===================

// get the means/variances/covariances in the respective buffers
// any one of the buffers can be NULL, in which case
// the corresponding moment will be skipped
// the length of buf_m, buf_v should be nvar->dim
// the length of buf_cov can be calculated using nvar_count_covariance()
// if ADD, stats numbers are added to the buffer numbers, weighted by WEIGHT,
// otherwise they replace whatever's there
void nvar_get_statistics(nvar_t *nvar, double *buf_m, double *buf_v, double *buf_cov,
			 int add, double weight)
{

    int i, j, idx = 0;

    for (i=0; i<nvar->dim; i++) {
	// calculate mean??
	if (buf_m) 
	    buf_m[i] = (add ?
			buf_m[i] + weight * nvar_calc_mean(nvar, i) :
			nvar_calc_mean(nvar, i));
	// calulate variance
	if (buf_v) 
	    buf_v[i] = (add ?
			buf_v[i] + weight * nvar_calc_covariance(nvar, i, i, FALSE) :
			nvar_calc_covariance(nvar, i, i, FALSE));
	// calculate covariances
	if (buf_cov) 
	    for (j=i+1; j<nvar->dim; j++, idx++) 
		buf_cov[idx] = (add ?
				buf_cov[idx] +
				weight * nvar_calc_covariance(nvar, i, j, FALSE) :
				nvar_calc_covariance(nvar, i, j, FALSE));
    }
    
}





// calculate the length of the vector needed to store the covariances
int nvar_count_covariances(nvar_t *nvar) {

    return nvar->dim * (nvar->dim - 1) / 2;
    
}
