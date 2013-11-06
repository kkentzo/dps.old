/*
=======
utils.h
=======

A suite of useful functions
 */

#include <glib.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#define MYMAX(A,B) ((A)>(B) ? (A) : (B))
#define MYMIN(A,B) ((A)<(B) ? (A) : (B))


#define LETTERS "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
#define LETTERS_SIZE 26

// MACRO to 
#define BFILL(BUF,LEN,V) {int _i; for (_i=0;_i<LEN;BUF[_i++]=V);}



//==========================================================================
//                     RANGE FUNCTIONS
//==========================================================================


//=======================================

typedef struct {

    double start;
    double stop;
    double step;

    double _cnt;

} range_t;



// creates and returns a range object
range_t *range_init(double start, double stop, double step);

// creates and returns a range *object by parsing a string
// with format "%f:%f:%f"
range_t *range_init_from_string(const char *st);


// free the range object
void range_free(range_t *range);


// copy the next value in *VAL - if no more values return FALSE
int range_next(range_t *range, double *val);


// reset the range
void range_reset(range_t *range);


// get number of elements in the range
int range_count(range_t *range);







//==========================================================================
//                     SORTING STUFF
//==========================================================================


// === compare functions ===

// for double elements
int dcmp_desc(const void *e1, const void *e2); 
int dcmp_asc(const void *e1, const void *e2); 
// for int elements
int icmp_desc(const void *e1, const void *e2); 
int icmp_asc(const void *e1, const void *e2); 


// sort the ARRAY (of length N - each element has size SIZE)
// the CMP function applies to the elements of ARRAY
// the function sorts the array's indices as well,
// which are stored in array IND
// !! IND must be initialized with the indices of ARRAY
void qsort_w_indices(void *array, int *ind, size_t n, size_t size,
		     int (*cmp)(const void *, const void *));



//==========================================================================
//                     NUMERICAL STUFF
//==========================================================================

// returns a random positive integer from /dev/urandom
int rand();

// returns a random unsigned integer from /dev/urandom
unsigned int urand();

// returns a random unsigned long from /dev/urandom
unsigned long lurand();

// define default values for mu and lamda in sigmoidal functions
#define SIG_MU 1.
#define SIG_LAMDA 1.

// sigmoid (logistic) function -- result in (0, lamda)
// lamda is the scaling factor, mu the steepness of the curve
double sigmoid0(double x, double mu, double lamda);

// sigmoid (logistic) function -- result in (-lamda/2, lamda/2)
// lamda is the scaling factor, mu the steepness of the curve
double sigmoid1(double x, double mu, double lamda);

// linear interpolation function
// returns the value of y(x) by interpolating between (x1,y1) and (x2,y2)
double interpolate(double x, double x1, double x2, double y1, double y2);

// convert the number of seconds to a nicely formatted string
// which is appended in st
void sec_to_str(long sec, GString *st);




//==========================================================================
//                            GSL STUFF
//==========================================================================

// returns the max/min/sum/mean of each column of the matrix in the returned gsl_vector
// similar to numpy's min(axis=0)
gsl_vector *gsl_matrix_max_across_rows(gsl_matrix *mat);
gsl_vector *gsl_matrix_min_across_rows(gsl_matrix *mat);
gsl_vector *gsl_matrix_sum_across_rows(gsl_matrix *mat);
gsl_vector *gsl_matrix_mean_across_rows(gsl_matrix *mat);


// returns the max/min/sum/mean of each row of the matrix in the returned gsl_vector
// similar to numpy's min(axis=1)
gsl_vector *gsl_matrix_max_across_cols(gsl_matrix *mat);
gsl_vector *gsl_matrix_min_across_cols(gsl_matrix *mat);
gsl_vector *gsl_matrix_sum_across_cols(gsl_matrix *mat);
gsl_vector *gsl_matrix_mean_across_cols(gsl_matrix *mat);



// add a value at an array element
void gsl_matrix_add_value_at_pos(gsl_matrix *mat, double val, int row, int col);
void gsl_vector_add_value_at_pos(gsl_vector *vec, double val, int pos);



//==========================================================================
//                         GPtrArray STUFF
//==========================================================================

// read a file and return the lines as an array of strings (contain the \n's)
GPtrArray *g_ptr_array_new_from_file(const char *fname);

// free all strings contained in arr and free arr as well
void g_ptr_array_free_with_func(GPtrArray *arr, GDestroyNotify func);



//==========================================================================
//                         GArray STUFF
//==========================================================================

// return an array with the specified range of integers in [lo, hi)
GArray *g_array_new_from_range(int lo, int hi);

// a destructor function for a GArray
void g_array_free_func(gpointer *data);



//==========================================================================
//                   A SYMMETRIC MATRIX CLASS
//                     e.g. a distance matrix
//==========================================================================
typedef struct {

    double **data;
    int size;

} sm_t;




// create a new sym. matrix with SIZE rows and SIZE columns
// and initialize all elements to the value FILL
sm_t *sm_new(int size, double fill);

// destroy the sym. matrix
void sm_free(sm_t *sm);

// fill the matrix with the supplied value
void sm_fill(sm_t *sm, double val);


// set an element of the matrix to a value
void sm_set(sm_t *sm, int row, int col, double val);

// sets VAL to be the value of the requested element
double sm_get(sm_t *sm, int row, int col);

// get the number of values in the symmetric matrix
// DIAGONAL determines whether to include or exclude the entries of the diagonal
int sm_count(sm_t *sm, gboolean diagonal);

// fill the supplied buffer with the values of the symmetric matrix
// CAREFUL: the length of the buffer should be equal to sm_count(sm, diagonal)
//         NO CHECK IS PERFORMED!!
// DIAGONAL determines whether to include or exclude the entries of the diagonal
void sm_flatten(sm_t *sm, double *buf, gboolean diagonal);

// fill the BUF with the entries of the diagonal
void sm_diagonal(sm_t *sm, double *buf);

// calculate and return the mean of the symmetric matrix
double sm_mean(sm_t *sm, gboolean diagonal);

// calculate and return the variance of the symmetric matrix
// the mean MU should also be supplied
double sm_variance(sm_t *sm, double mu, gboolean diagonal);

// print the matrix 
void sm_print(sm_t *sm);


//==========================================================================
//                           STRING STUFF
//==========================================================================


// creates a gstring with all the file's contents 
// supports gzipped files as well!
GString *g_string_new_from_file(const char *fname_);

// creates a gstring by copying the contents of str in the specified range
GString *g_string_new_from_str_with_range(const char *init, int from, int to);

// creates a gstring by joining other gstrings
GString *g_string_new_by_joining(GPtrArray *arr, const char sep,
				 gboolean begin_with_sep);

// creates a new string from VAL that corresponds to a letter of the alphabet
// if VAL is > 26, then a number is appended
GString *g_string_new_itol(int val);

// a destructor function for a gstring
void g_string_free_func(gpointer st);

// remove any leading or trailing whitespace or new line character from st
GString *g_string_strip(GString *st);

// remove any trailing whitespace or new line character from st
GString *g_string_rstrip(GString *st);

// remove any leading whitespace or new line character from st
GString *g_string_lstrip(GString *st);

// split the string using the character separator and return an array of strings
// remember to free the GPtrArray using : 
// g_ptr_array_free_with_func(ARRAY_NAME, g_string_free_func);
GPtrArray *g_string_split(GString *st, const char sep);

// tokenize a string (C version)
// returns an array of GString pointers (tokens) by parsing the nested list
// eg tokenize('((a b) (c d) foo)')
// yields {'(a b)', '(c d)', 'foo'}
// remember to free the GPtrArray using : 
// g_ptr_array_free_with_func(ARRAY_NAME, g_string_free_func);
GPtrArray *g_string_tokenize(GString *st);

// add the path to gstring 
GString *g_string_append_path(GString *st, const char *path);


//==========================================================================
//                        FILE-RELATED STUFF
//==========================================================================

// compress the file FNAME
// will produce the file FNAME.gz
int compress_file(const char *fname);

// decompress the file FNAME in the same directory as FNAME
int decompress_file(const char *fname);


int compress_dir(const char *fname, const char fmt, gboolean delete);

// check if a file exists
gboolean fexists(const char *path);









// THE SELECTION STUFF BELOW IS NOW OBSOLETE - USE sampling.h

//==========================================================================
//                     SELECTION STUFF 
//==========================================================================


// // select randomly n unique integers in [lo, hi) and store them in selected
// void uselect(int *selected, int n, int lo, int hi, void *rng);


// // ========== roulette wheel selection ==========
// // returns an index of vec based on the probabilities stored in vec
// // CAREFUL :: sum(probs) should be 1 - no check is performed
// int roulette_p(const double *vec, int len, void *rng);


// // returns an index of vec based on the values (not probs) stored in vec
// // calculates the sum of vec
// int roulette_v(const double *vec, int len, void *rng);

// // returns an index of vec based on the values (not probs) stored in vec
// // provide also the sum of vec
// int roulette_v_s(const double *vec, int len, double vec_sum, void *rng);

// // select an integer from IVEC, based on the values stored in VEC
// // the integers in IVEC index VEC, so the length of VEC must be equal to
// // the largest possible integer stored in IVEC
// // however only the entries of VEC that correspond to indices in IVEC must
// // hold meaningful values
// // VEC_SUM is the sum of VEC (used to calculate selection probabilities)
// // if REMOVE, then the selected index is removed from IVEC
// int select_index(GArray *ivec, double *vec, double vec_sum, gboolean remove, void *rng);



// // ===============================================================================
// // the two functions below are based on an algorithm found in:
// // http://stackoverflow.com/questions/6443176/how-can-i-generate-a-random-number-within-a-range-but-exclude-some
// // ===============================================================================

// // return a random int in [lo, hi) excluding the value EVAL
// int get_rnd_int_excluding_value(int lo, int hi, int eval, void *rng);

// // return a random int in [lo, hi) excluding the values in array EVALS
// // ** elements in EVALS should be in ascending order and mutually exclusive **
// int get_rnd_int_excluding_values(int lo, int hi, int *evals, int evals_len, void *rng);





// ================================================================
// return a double from a uniform distribution in [lo, hi]
//double get_rnd_double_with_bounds(double lo, double hi, void *rng);
