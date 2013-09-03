/*
sampling.h

Various sampling functions 

** sampling from an integer sequence with/without replacement
** sampling from an integer sequence by excluding specific items

 */


//===============================================================================
//                       UTILITY FUNCTIONS
//===============================================================================


// mutate the supplied value by sampling from a uniform distribution in [VAL-DELTA, VAL+DELTA]
// makes sure that the returned value is within the supplied LO and HI and boundaries
// PBC indicates whether there exist periodic boundary conditions
double mutate(double val, double delta, double lo, double hi, int pbc, void *rng);


// normalizes the vector V (of length N) so that \sum_i v_i = 1
// if VSUM > 0, it is used as the vector's sum 
// otherwise the sum of the vector is computed
void vnormalize(double *v, int n, double vsum);



// shuffles (randomizes the order of) the contents of the VEC
void shuffle(void *v, int n, int size, void *rng);




//=================================
//     SAMPLING - NO EXCLUSION
//=================================

// sample k integers from [0,n) and store them in dest
// if src==NULL, all options are considered with equal probability
// otherwise, each option i is considered with probability src[i]
// As such it must be that len(dest) == k and len(src) == n
// REPLACE determines whether sampling occurs with replacement or not
// NORM determines whether src needs to be normalized (so that \sum src_i = 1)
void sample(int k, int n, int *dest, double *src, int replace, int norm, void *rng);



//=================================
//   SAMPLING - WITH EXCLUSION
//=================================


// sample k integers from [0,n) and store them in dest
// sampling excludes the values stored in ev (with length en)
// all options are considered with equal probability
// based on an algorithm found at :
// http://stackoverflow.com/questions/6443176/how-can-i-generate-a-random-number-within-a-range-but-exclude-some
void sample_excluding(int k, int n, int *dest, int *ev, int en, void *rng);





