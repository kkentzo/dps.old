/*
nvar.h

n-variate descriptive statistics

*/





typedef struct {

    // the number of variables
    int dim;

    // the number of supplied observations (sum of the weights)
    double n;

    // vector with the weighted sums of each variables (length: dim)
    double *wsums; 

    // symmetric matrix (sm_t *) with the sums of the weighted pair-wise products
    void *sum_wprods;

    // variable labels
    char **labels;

    // pairwise variable labels (for covariance labels)
    // array size can be calculated using nvar_count_covariances()
    // it's actually == n * (n-1) / 2
    char **plabels;

} nvar_t;




nvar_t *nvar_init(int dim);

void nvar_free(nvar_t *nvar);

void nvar_reset(nvar_t *nvar);


// add variable labels (copy the strings)
void nvar_assign_labels(nvar_t *nvar, const char **labels);

// print the labels
void nvar_print_cov_labels(nvar_t *nvar);

// add the data point specified in VALUES (of length nvar->dim) 
// with the supplied WEIGHT
void nvar_add(nvar_t *nvar, double *values, double weight);


// add the supplied stats (NVAR2) to NVAR1
void nvar_add_stats(nvar_t *nvar1, nvar_t *nvar2);


// remove the data point specified in VALUES (of length nvar->dim) 
// with the supplied WEIGHT
void nvar_remove(nvar_t *nvar, double *values, double weight);


// ================= CALCULATE STATISTICS ===================

// calculate the mean of the specified variable
double nvar_calc_mean(nvar_t *nvar, int idx);


// calculate the covariance between variables IDX1 and IDX2
// if IDX1 == IDX2 returns the variable's variance 
// if (NORM) then returns the correlation coefficient between IDX1, IDX2
//            and if IDX1==IDX2 then the result is 1
double nvar_calc_covariance(nvar_t *nvar, int idx1, int idx2, int norm);


// get the means/variances/covariances in the respective buffers
// any one of the buffers can be NULL, in which case
// the corresponding moment will be skipped
// the length of buf_m, buf_v should be nvar->dim
// the length of buf_cov can be calculated using nvar_count_covariance()
// if ADD, stats numbers are added to the buffer numbers, weighted by WEIGHT,
// otherwise they replace whatever's there
void nvar_get_statistics(nvar_t *nvar, double *buf_m, double *buf_v, double *buf_cov,
			 int add, double weight);



// calculate the length of the vector needed to store the covariances
int nvar_count_covariances(nvar_t *nvar);




// ================= HDF5 FUNCTIONS ======================

