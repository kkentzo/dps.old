#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <unistd.h> // for chdir()

#include <glib.h>
#include <gsl/gsl_rng.h>

#include "clib.h"



//==========================================================================
//                            RANGE STUFF
//==========================================================================

// creates and returns a range object
range_t *range_init(double start, double stop, double step) {

    range_t *r = malloc(sizeof(range_t));

    r->start = start;
    r->stop = stop;
    r->step = step;

    r->_cnt = 0;

    return r;

}


// creates and returns a range *object by parsing a string
// with format "%f:%f:%f"
range_t *range_init_from_string(const char *st) {

    // create a gstring from st
    GString *gst = g_string_new(st);

    // split the gstring
    GPtrArray *tokens = g_string_split(gst, ':');

    assert(tokens->len == 3);

    // parse the values

    GString *s = g_ptr_array_index(tokens, 0);
    double start = atof(s->str);
    g_string_free(s, TRUE);

    s = g_ptr_array_index(tokens, 1);
    double stop = atof(s->str);
    g_string_free(s, TRUE);

    s = g_ptr_array_index(tokens, 2);
    double step = atof(s->str);
    g_string_free(s, TRUE);

    // free the pointer array (strings have been freed)
    g_ptr_array_free(tokens, TRUE);

    // free the gstring
    g_string_free(gst, TRUE);

    return range_init(start, stop, step);

}






// free the range object
void range_free(range_t *range) {

    if (range)
	free(range);

}




// copy the next value in *VAL - if no more values return FALSE
int range_next(range_t *range, double *val) {

    // calculate value
    double v = range->start + range->step * range->_cnt;
    if (v >= range->stop)
	return FALSE;
    else {
	range->_cnt ++;
	*val = v;
	return TRUE;
    }
	
    
}



// reset the range
void range_reset(range_t *range) {

    range->_cnt = 0;
    
}




// get number of elements in the range
int range_count(range_t *range) {

    return (range->stop - range->start) / range->step;

}





//==========================================================================
//                     SORTING STUFF
//==========================================================================


// === compare functions ===

// for double elements - ascending order
int dcmp_asc(const void *e1, const void *e2) {

    double v1 = *(double *)e1;
    double v2 = *(double *)e2;

    return (v1 - v2 > 0 ? 1 : (v1 - v2 < 0 ? -1 : 0));
    
}




// for double elements - descending order
int dcmp_desc(const void *e1, const void *e2) {
    

    double v1 = *(double *)e1;
    double v2 = *(double *)e2;

    return (v1 - v2 > 0 ? -1 : (v1 - v2 < 0 ? 1 : 0));
    
}



// for int elements - ascending order
int icmp_asc(const void *e1, const void *e2) {

    int v1 = *(int *)e1;
    int v2 = *(int *)e2;

    return v1 - v2;

}



// for int elements - descending order
int icmp_desc(const void *e1, const void *e2) {

    int v1 = *(int *)e1;
    int v2 = *(int *)e2;

    return v2 - v1;
    
}



// sort the ARRAY (of length N - each element has size SIZE)
// the CMP function applies to the elements of ARRAY
// the function sorts the array's indices as well,
// which are stored in array IND
// !! IND must be initialized with the indices of ARRAY
void qsort_w_indices(void *array, int *ind, size_t n, size_t size,
		     int (*cmp)(const void *, const void *)) {

    int _cmp_(const void *e1, const void *e2) {

	int i1 = *(int *)e1;
	int i2 = *(int *)e2;

	return cmp((const void *)(array + i1 * size),
		   (const void *)(array + i2 * size));
	
    }

    // sort the ind array
    qsort(ind, n, sizeof(int), _cmp_);

    // copy the array in buf
    void *buf = malloc(n * size);
    memmove(buf, array, n * size);
    // rearrange array items acc. to sorted indices
    int i;
    for (i=0; i<n; i++)
	memcpy(array+i*size, buf+ind[i]*size, size);
    
    free(buf);
    
}






//==========================================================================
//                     NUMERICAL STUFF
//==========================================================================





// returns a random positive integer from /dev/urandom
int rand() {

    int val;
    FILE *fp = fopen("/dev/urandom", "rb");
    fread(&val, sizeof val, 1, fp);
    fclose(fp);

    return abs(val);
    
}





// returns a random unsigned integer from /dev/urandom
unsigned int urand() {

    unsigned int val;
    FILE *fp = fopen("/dev/urandom", "rb");
    fread(&val, sizeof val, 1, fp);
    fclose(fp);

    return val;
    
}





// returns a random unsigned long from /dev/urandom
unsigned long lurand() {

    unsigned long val;
    FILE *fp = fopen("/dev/urandom", "rb");
    fread(&val, sizeof val, 1, fp);
    fclose(fp);

    return val;
    
}    







// sigmoid (logistic) function -- result in (0, lamda)
// lamda is the scaling factor, mu the steepness of the curve
double sigmoid0(double x, double mu, double lamda) {

    return lamda / (1 + exp(- mu * x));
    
}






// sigmoid (logistic) function -- result in (-lamda/2, lamda/2)
// lamda is the scaling factor, mu the steepness of the curve
double sigmoid1(double x, double mu, double lamda) {

    double h = exp(-mu*x);

    return lamda * (1 - h) / (2 * (1 + h));

    
}



// linear interpolation function
// returns the value of y(x) by interpolating between (x1,y1) and (x2,y2)
double interpolate(double x, double x1, double x2, double y1, double y2) {

    return y1 + (y2-y1) * (x-x1) / (x2-x1);
    
}





// convert the number of seconds to a nicely formatted string
// which is appended in st
void sec_to_str(long sec, GString *st) {

    long days, hours, mins;
    ldiv_t res;

    sec = labs(sec);

    //NSLog(@"%ld", sec);
    
    res = ldiv(sec, 60);
    mins = res.quot;
    sec = res.rem;

    res = ldiv(mins, 60);
    hours = res.quot;
    mins = res.rem;

    res = ldiv(hours, 24);
    days = res.quot;
    hours = res.rem;

    g_string_append_printf(st, "%ldd %02ld:%02ld:%02ld", days, hours, mins, sec);
    
}





//==========================================================================
//                            GSL STUFF
//==========================================================================

// returns the max/min of each column of the matrix in the returned gsl_vector
// similar to numpy's min(axis=0)
gsl_vector *gsl_matrix_max_across_rows(gsl_matrix *mat) {

    int i;
    gsl_vector *vec = gsl_vector_alloc(mat->size2);
    gsl_vector_view column;

    for (i=0; i<vec->size; i++) {
	// get view of i^th column
	column = gsl_matrix_column(mat, i);
	// find max of column and add it to vec
	gsl_vector_set(vec, gsl_vector_max(&column.vector), i);
    }

    return vec;
    
}



gsl_vector *gsl_matrix_min_across_rows(gsl_matrix *mat) {

    int i;
    gsl_vector *vec = gsl_vector_alloc(mat->size2);
    gsl_vector_view column;

    for (i=0; i<vec->size; i++) {
	// get view of i^th column
	column = gsl_matrix_column(mat, i);
	// find min of column and add it to vec
	gsl_vector_set(vec, gsl_vector_min(&column.vector), i);
    }

    return vec;
    
}





gsl_vector *gsl_matrix_sum_across_rows(gsl_matrix *mat) {

    int row, col;
    gsl_vector *vec = gsl_vector_alloc(mat->size2);
    gsl_vector_set_zero(vec);


    for (col=0; col<mat->size2; col++)
	for (row=0; row<mat->size1; row++)
	    gsl_vector_add_value_at_pos(vec, gsl_matrix_get(mat, row, col), col);

    return vec;
    
}





gsl_vector *gsl_matrix_mean_across_rows(gsl_matrix *mat) {

    int row, col;
    gsl_vector *vec = gsl_vector_alloc(mat->size2);
    gsl_vector_set_zero(vec);

    for (col=0; col<mat->size2; col++) {    
	for (row=0; row<mat->size1; row++) 
	    gsl_vector_add_value_at_pos(vec, gsl_matrix_get(mat, row, col), col);

	// divide to get the mean
	gsl_vector_set(vec, gsl_vector_get(vec, col) / mat->size1, col);
	
    }

    return vec;
    
}







// returns the max/min of each row of the matrix in the returned gsl_vector
// similar to numpy's min(axis=1)
gsl_vector *gsl_matrix_max_across_cols(gsl_matrix *mat) {

    int i;
    gsl_vector *vec = gsl_vector_alloc(mat->size1);
    gsl_vector_view row;

    for (i=0; i<vec->size; i++) {
	// get view of i^th row
	row = gsl_matrix_row(mat, i);
	// find max of row and add it to vec
	gsl_vector_set(vec, gsl_vector_max(&row.vector), i);
    }

    return vec;
    
}



gsl_vector *gsl_matrix_min_across_cols(gsl_matrix *mat) {

    int i;
    gsl_vector *vec = gsl_vector_alloc(mat->size1);
    gsl_vector_view row;

    for (i=0; i<vec->size; i++) {
	// get view of i^th row
	row = gsl_matrix_row(mat, i);
	// find max of row and add it to vec
	gsl_vector_set(vec, gsl_vector_min(&row.vector), i);
    }

    return vec;
    
}






gsl_vector *gsl_matrix_sum_across_cols(gsl_matrix *mat) {

    int row, col;
    gsl_vector *vec = gsl_vector_alloc(mat->size1);
    gsl_vector_set_zero(vec);


    for (row=0; row<mat->size1; row++)    
	for (col=0; col<mat->size2; col++)
	    gsl_vector_add_value_at_pos(vec, gsl_matrix_get(mat, row, col), row);

    return vec;
    
}






gsl_vector *gsl_matrix_mean_across_cols(gsl_matrix *mat) {

    int row, col;
    gsl_vector *vec = gsl_vector_alloc(mat->size1);
    gsl_vector_set_zero(vec);

    for (row=0; row<mat->size1; row++) {
	for (col=0; col<mat->size2; col++)
	    gsl_vector_add_value_at_pos(vec, gsl_matrix_get(mat, row, col), row);

	// divide to get the mean
	gsl_vector_set(vec, gsl_vector_get(vec, row) / mat->size2, row);
	
    }

    return vec;
    
}




// add a value at an array element
void gsl_matrix_add_value_at_pos(gsl_matrix *mat, double val, int row, int col) {

    gsl_matrix_set(mat, row, col, gsl_matrix_get(mat, row, col) + val);
    
}



void gsl_vector_add_value_at_pos(gsl_vector *vec, double val, int pos) {

    gsl_vector_set(vec, pos, gsl_vector_get(vec, pos) + val);

}



//==========================================================================
//                     GPtrArray STUFF
//==========================================================================

// read a file and return the lines as an array of strings
GPtrArray *g_ptr_array_new_from_file(const char *fname) {

    FILE *f = fopen(fname, "r");
    if (! f)
	return NULL;

    //GPtrArray *arr = g_ptr_array_new_with_free_func(g_string_free_func);
    GPtrArray *arr = g_ptr_array_new();

    char *line = NULL;
    size_t len = 0;
    ssize_t bytes_read;

    while ( ( bytes_read = getline(&line, &len, f)) != -1 )
	// if the line is not empty
	if (bytes_read > 0) 
	    // allocate space for string and add it to pointer array
	    g_ptr_array_add(arr, g_string_new_len(line, len));

    // free the line pointer
    free(line);
    // close the file
    fclose(f);

    return arr;
    
}





// free all objects contained in arr and free arr as well
void g_ptr_array_free_with_func(GPtrArray *arr, GDestroyNotify func) {

    int i;
    // free all objects stored in arr
    for (i=0; i<arr->len; i++)
	func(g_ptr_array_index(arr, i));

    // free arr
    g_ptr_array_free(arr, TRUE);
	     
}








//==========================================================================
//                         GArray STUFF
//==========================================================================

// return an array with the specified range of integers in [lo, hi)
GArray *g_array_new_from_range(int lo, int hi) {

    assert(lo < hi);

    GArray *garr = g_array_sized_new(FALSE, FALSE, sizeof(int), hi-lo);
    int i;

    for (i=lo; i<hi; i++)
	g_array_append_val(garr, i);

    return garr;
    
}



// a destructor function for a GArray
void g_array_free_func(gpointer *data) {

    g_array_free((GArray *)data, TRUE);

}




//==========================================================================
//                   A SYMMETRIC SQUARE MATRIX CLASS
//==========================================================================

// create a new sym. matrix with LEN rows and LEN columns
// and initialize all elements to the value FILL
sm_t *sm_new(int size, double fill) {

    // allocate space for the structure
    sm_t *sm = malloc(sizeof(sm_t));
    // set size
    sm->size = size;

    // allocate space for the array of pointers
    // to the rows of the matrix
    sm->data = malloc(sm->size * sizeof(double *));

    // allocate space for each row of the array
    int i, j;
    for (i=0; i<sm->size; i++) {
	sm->data[i] = malloc((sm->size - i) * sizeof(double));
	for (j=0; j<sm->size-i; j++)
	    sm->data[i][j] = fill;
    }


    return sm;
    
}

// destroy the sym. matrix
void sm_free(sm_t *sm) {

    if (! sm)
	return;

    // free all rows
    int i;
    for (i=0; i<sm->size; i++)
	free(sm->data[i]);
    // free the array of row pointers
    free(sm->data);
    // free the struct object
    free(sm);
    
}




// fill the matrix with the supplied value
void sm_fill(sm_t *sm, double val) {

    int i,j;

    for (i=0; i<sm->size; i++)
	for (j=0; j<sm->size-i; j++)
	    sm->data[i][j] = val;
    
}



// set an element of the matrix to a value
void sm_set(sm_t *sm, int row, int col, double val) {

    // swap row and column??
    if (row > col) {
	int tmp = row; row = col; col = tmp;
    }

    // set value
    sm->data[row][col-row] = val;
    
}



// sets VAL to be the value of the requested element
double sm_get(sm_t *sm, int row, int col) {

    // swap row and column??
    if (row > col) {
	int tmp = row; row = col; col = tmp;
    }

    // retrieve and return value
    return sm->data[row][col-row];
    
}




// get the number of values in the symmetric matrix
int sm_count(sm_t *sm, gboolean diagonal) {

    return (diagonal ? sm->size : 0) + sm->size * (sm->size-1) / 2;
    
}





// fill the supplied buffer with the values of the symmetric matrix
// CAREFUL: the length of the buffer should be equal to sm_count(sm)
//         NO CHECK IS PERFORMED!!
// DIAGONAL determines whether to include or exclude the entries of the diagonal
void sm_flatten(sm_t *sm, double *buf, gboolean diagonal) {

    int i,j;
    int cnt = 0;

    for (i=0; i<sm->size; i++)
	for (j=(diagonal ? i : i+1); j<sm->size; j++) {
	    buf[cnt] = sm_get(sm, i, j);
	    cnt += 1;
	}

}



// fill the BUF with the entries of the diagonal
void sm_diagonal(sm_t *sm, double *buf) {

    int i;
    for (i=0; i<sm->size; i++)
	buf[i] = sm_get(sm, i, i);
    
}





// calculate and return the mean of the symmetric matrix
double sm_mean(sm_t *sm, gboolean diagonal) {

    int i,j;
    double sum = 0;
    int cnt = 0;

    for (i=0; i<sm->size; i++)
	for (j=(diagonal ? i : i+1); j<sm->size; j++) {
	    sum += sm_get(sm, i, j);
	    cnt += 1;
	}

    // the statement below works
    //assert(cnt == sm_count(sm));

    return sum / cnt;
    
}





// calculate and return the variance of the symmetric matrix
double sm_variance(sm_t *sm, double mu, gboolean diagonal) {

    int i,j;
    double sum = 0;
    double x;
    int cnt = 0;

    for (i=0; i<sm->size; i++)
	for (j=(diagonal ? i : i+1); j<sm->size; j++) {
	    x = sm_get(sm, i, j);
	    sum += (x - mu) * (x - mu);
	    cnt += 1;
	}

    return sum / (cnt - 1);
    
}




// print the matrix 
void sm_print(sm_t *sm) {

    int i, j;
    for (i=0; i<sm->size; i++) {
	for (j=0; j<sm->size; j++)
	    printf("%5.1f ", sm_get(sm, i, j));
	printf("\n");
    }
    
}






//==========================================================================
//                     STRING STUFF
//==========================================================================



// returns a gstring with all the file's contents
// supports gzipped files as well!
GString *g_string_new_from_file(const char *fname_) {

    gboolean gz = g_str_has_suffix(fname_, ".gz");

    // from now on, we're working with the GString file name
    GString *fname = g_string_new(fname_);

    if (gz) {
	// decompress the file (deletes the gzipped file)
	decompress_file(fname->str);
	// remove the .gz extension from the file name
	g_string_truncate(fname, fname->len - 3);
    }

    // load and read file into string
    FILE *fp = fopen(fname->str, "r");
    // did file open succeed??
    if (! fp) 
	return NULL;

    // create string of file contents
    GString *contents = g_string_new(NULL);

    // read file into contents string
    char ch;
    while ((ch = fgetc(fp)) != EOF)
	g_string_append_c(contents, ch);

    // close file
    fclose(fp);

    if (gz)
    	// recompress file
    	compress_file(fname->str);

    g_string_free(fname, TRUE);

    return contents;

}    




// creates a gstring by copying the contents of str in the specified range
GString *g_string_new_from_str_with_range(const char *init, int from, int to) {

    if (from >= to)
	return NULL;

    int len_init = strlen(init);
    // make sure we don't overstep init's length
    if (to >= len_init) 
	to = len_init;

    GString *st = g_string_sized_new(len_init);
    int i;
    for (i=from; i<to; i++)
	g_string_append_c(st, init[i]);

    return st;
    
}





// creates a gstring by joining other gstrings
GString *g_string_new_by_joining(GPtrArray *arr, const char sep,
				 gboolean begin_with_sep) {

    if (arr == NULL || arr->len == 0)
	return NULL;

    GString *s, *st = g_string_new(NULL);

    if (begin_with_sep)
	g_string_append_c(st, sep);
    
    int i;
    for (i=0; i<arr->len; i++) {
	s = g_ptr_array_index(arr, i);
	g_string_append(st, s->str);
	if (i < arr->len - 1)
	    g_string_append_c(st, sep);
	    
    }

    return st;
    
}




// creates a new string from VAL that corresponds to a letter of the alphabet
// if VAL is > 26, then a number is appended
GString *g_string_new_itol(int val) {

    assert(val >= 0);

    GString *st = g_string_new("");

    // write contents
    if (val >= LETTERS_SIZE)
	g_string_printf(st, "%c%d", LETTERS[val % LETTERS_SIZE], val / LETTERS_SIZE);
    else
	g_string_printf(st, "%c", LETTERS[val]);

    return st;
}




// a destructor function for a gstring 
void g_string_free_func(gpointer st) {

    g_string_free((GString *) st, TRUE);

}






// remove any leading or trailing whitespace or new line character from st
GString *g_string_strip(GString *st) {

    g_string_lstrip(st);
    g_string_rstrip(st);
    return st;
    
}




// remove any trailing whitespace or new line character from st
GString *g_string_rstrip(GString *st) {

    int pos = st->len;
    
    while (pos--)
	// isspace includes space, formfeed, newline, cr, tab, vtab
	if (! isspace(st->str[pos]))
	    break;

    return g_string_set_size(st, pos+1);

}





// remove any leading whitespace or new line character from st
GString *g_string_lstrip(GString *st) {

    int i;

    for (i=0; i<st->len; i++)
	// isspace includes space, formfeed, newline, cr, tab, vtab
	if (! isspace(st->str[i]))
	    break;

    // erase the first i characters
    g_string_erase(st, 0, i);

    return st;
}






// split the string using the separator and return an array of strings
// remember to free the GPtrArray using : 
// g_ptr_array_free_with_func(ARRAY_NAME, g_string_free_func);
GPtrArray *g_string_split(GString *st, const char sep) {

    int i, marker = 0;
    GPtrArray *arr = g_ptr_array_new();
    GString *s;
    
    for (i=0; i<st->len; i++)
    	if (st->str[i] == sep) {
	    // got it! - create string
	    s = g_string_new_from_str_with_range(st->str,
						 marker,
						 i);
	    // set marker to the next character - if the next char
	    // happens to conicide with sep, then s will be NULL
	    marker = i+1;
	    if (s) 
		g_ptr_array_add(arr, s);

	}

    // add any remainder to the pointer array
    if (marker < st->len)
	g_ptr_array_add(arr, g_string_new_from_str_with_range(st->str,
							      marker,
							      st->len));

    return arr;
    
}



// tokenize a string (C version)
// returns an array of GString pointers (tokens) by parsing the nested list
// eg tokenize('((a b) (c d) foo)')
// yields {'(a b)', '(c d)', 'foo'}
// remember to free the GPtrArray using : 
// g_ptr_array_free_with_func(ARRAY_NAME, g_string_free_func);
GPtrArray *g_string_tokenize(GString *st) {

    GString *get_next_token(int *pos) {
	
	int ctr=0;
	char ch;
	// create token
	GString *token = g_string_new(NULL);

	while ( (*pos) < st->len - 1 ) {
	    ch = st->str[*pos];
	    if (ch == '(')
		ctr += 1;
	    if (ch == ')')
		ctr -= 1;
	    g_string_append_c(token, ch);
	    (*pos)++;
	    if (ctr == 0 && isspace(ch)) 
		return g_string_rstrip(token);
	}

	return g_string_rstrip(token);

    }

    // check that the supplied string begins and ends with a parenthesis
    if (st->len <= 2)
	return NULL;
    if ( ! (st->str[0] == '(' && st->str[st->len - 1] == ')') )
	return NULL;

    // initialize the array of tokens
    //GPtrArray *tokens = g_ptr_array_new_with_free_func(g_string_free_func);
    GPtrArray *tokens = g_ptr_array_new();
    // initialize the position
    int pos = 1; // ignore the opening parenthesis
    while (pos < st->len - 1) { // ignore the closing parenthesis
	g_ptr_array_add(tokens, get_next_token(&pos));
    }

    return tokens;

}





// add the path to gstring 
GString *g_string_append_path(GString *st, const char *path) {

    // if st is empty, just append path and return
    if (st->len == 0)
	return g_string_append(st, path);

    int len_path = strlen(path);
    // is path empty??
    if (len_path == 0)
	return st;


    // do we have a trailing and a leading slash joined?
    if (st->str[st->len-1] == '/' && path[0] == '/')
    // remove one '/'
	g_string_truncate(st, st->len-1);
    
    // do we have no slashes to join?
    else if (st->str[st->len-1] != '/' && path[0] != '/')
	// add one '/'
	g_string_append_c(st, '/');

    // append the path to string
    g_string_append(st, path);

    return st;
    
}




//==========================================================================
//                           VARIOUS STUFF
//==========================================================================

// replaces FNAME with the compressed version FNAME.gz
int compress_file(const char *fname) {

    if (! fexists(fname))
	return -1;

    // form compress command
    GString *cmd = g_string_new(NULL);
    g_string_printf(cmd, "gzip -f %s", fname);

    // and execute
    int status = system(cmd->str);

    // free cmd string
    g_string_free(cmd, TRUE);

    return status;

}




// decompress the file FNAME in the same directory as FNAME
int decompress_file(const char *fname) {

    if (! fexists(fname))
	return -1;

    // form decompress command
    GString *cmd = g_string_new(NULL);
    g_string_printf(cmd, "gzip -d %s", fname);

    // and execute
    int status = system(cmd->str);

    // free cmd string
    g_string_free(cmd, TRUE);

    return status;
    
}








// compress the file or directory PATH
// using the file format FMT ['j'|'z']
int compress_dir(const char *path, const char fmt, gboolean delete) {

    if (strlen(path) == 0 || path[0] == '~')
	return -1;

    int status;
    
    // form path as a gstring
    GString *fpath = g_string_new(path);

    // get path components
    GPtrArray *pcomponents = g_string_split(fpath, '/');
    // get name of file/directory
    GString *name = g_ptr_array_remove_index(pcomponents, pcomponents->len-1);
    // get parent directory
    GString *parent_dir = g_string_new_by_joining(pcomponents, '/',
						  fpath->str[0] == '/');
    if (parent_dir == NULL)
	parent_dir = g_string_new("./");

    // form compress command 
    GString *cmd = g_string_new(NULL);
    g_string_printf(cmd, "tar c%cvf %s.%s %s", fmt, name->str, 
		    (fmt=='j' ? "tar.bz2" : "tar.gz"), name->str);

    // == start performing actions ===

    // get current working directory
    char *cwd = getcwd(NULL, 0);
    printf("cwd=%s\n", cwd);

    // change directory to file's parent directory
    status = chdir(parent_dir->str);
    if (status)
	goto cleanup;

    // and execute
    status = system(cmd->str);

    if (status)
	goto cleanup;
    
    if (delete) {
	// delete file
	g_string_printf(cmd, "rm -rf %s", name->str);
	status = system(cmd->str);
    }


 cleanup:

    // restore cwd directory
    chdir(cwd);
    free(cwd);
    
    // free gstrings
    g_string_free(fpath, TRUE);
    g_string_free(name, TRUE);
    g_string_free(parent_dir, TRUE);
    g_string_free(cmd, TRUE);

    // free array of path components 
    g_ptr_array_free_with_func(pcomponents, g_string_free_func);

    return status;

}



// check if a file exists
gboolean fexists(const char *path) {

    FILE *fp = fopen(path, "r");
    gboolean file_exists = (fp != NULL);

    if (fp)
	fclose(fp);

    return file_exists;
    
}









//==========================================================================
//                           OBSOLETE STUFF
//==========================================================================



// returns the sum of the vector
/* double vector_sum(const double *vec, int len) { */

/*     double sum = 0; */
/*     int i; */
/*     for (i=0; i<len; i++) */
/* 	sum += vec[i]; */
/*     return sum; */
/* } */





/* // select randomly n unique integers in [lo, hi) and store them in selected */
/* void uselect(int *selected, int n, int lo, int hi, void *rng) { */

/*     // determine the length of the integer sequence */
/*     int len = hi - lo; */
/*     assert(len > 0 && n <= len); */

/*     int i, j; */

/*     // should we return all the numbers in the sequence?? */
/*     if (n == len) { */
/* 	for (i=0; i<n; i++) */
/* 	    selected[i] = lo + i; */
/* 	return; */
/*     } */

/*     // populate sequence of integers */
/*     GArray *sequence = g_array_new_from_range(lo, hi); */

/*     // pop n values from array and store them in SELECTED */
/*     for (i=0; i<n; i++) { */
/* 	// select a sequence index */
/* 	j = gsl_rng_uniform_int((gsl_rng *)rng, sequence->len); */
/* 	// store that number in selected */
/* 	selected[i] = g_array_index(sequence, int, j); */
/* 	// remove that number from the sequence */
/* 	g_array_remove_index_fast(sequence, j); */
/*     } */
	
/*     g_array_free(sequence, TRUE); */

/* } */









/* // ========== roulette wheel selection ========== */

/* // returns an index of vec based on the probabilities stored in vec */
/* // CAREFUL :: sum(probs) should be 1 - no check is performed */
/* int roulette_p(const double *vec, int len, void *rng) { */

/*     int i; */
/*     // draw a random number in [0,1] */
/*     double r = gsl_rng_uniform((gsl_rng *)rng); */
/*     // spin the wheel */
/*     double csum = 0.; */
/*     for (i=0; i < len; i++) { */
/* 	csum += vec[i]; */
/* 	if (r < csum) */
/* 	    return i; */
/*     } */

/*     // return last vector index */
/*     return --i; */
    
/* } */




/* // returns an index of vec based on the values (not probs) stored in vec) */
/* // calculates the sum of vec */
/* int roulette_v(const double *vec, int len, void *rng) { */

/*     // calculate sum of vec */
/*     double vsum = vector_sum(vec, len); */
/*     int i; */
/*     // draw a random number in [0,1] */
/*     double r = gsl_rng_uniform((gsl_rng *)rng); */
/*     // spin the wheel */
/*     double csum = 0.; */
/*     for (i=0; i < len; i++) { */
/* 	// calculate probability of i^th element */
/* 	csum += vec[i] / vsum; */
/* 	if (r < csum) */
/* 	    return i; */
/*     } */

/*     // return last vector index */
/*     return --i; */
    
/* } */




/* // returns an index of vec based on the values (not probs) stored in vec */
/* // provide also the sum of vec */
/* int roulette_v_s(const double *vec, int len, double vec_sum, void *rng) { */

/*     int i; */
/*     // draw a random number in [0,1] */
/*     double r = gsl_rng_uniform((gsl_rng *)rng); */
/*     // spin the wheel */
/*     double csum = 0.; */
/*     for (i=0; i < len; i++) { */
/* 	// calculate probability of i^th element */
/* 	csum += vec[i] / vec_sum; */
/* 	if (r < csum) */
/* 	    return i; */
/*     } */

/*     // return last vector index */
/*     return --i; */
    
/* } */





/* // select an integer from IVEC, based on the values stored in VEC */
/* // the integers in IVEC index VEC, so the length of VEC must be equal to */
/* // the largest possible integer stored in IVEC */
/* // however only the entries of VEC that correspond to indices in IVEC must */
/* // hold meaningful values */
/* // VEC_SUM is the sum of VEC (used to calculate selection probabilities) */
/* // if REMOVE, then the selected index is removed from IVEC */
/* int select_index(GArray *ivec, double *vec, double vec_sum, gboolean remove, void *rng) { */

/*     assert(ivec->len > 0); */

/*     int i; // indexes IVEC */
/*     int idx; // indexes VEC (values contained in IVEC) */

/*     // draw a random number in [0,1] */
/*     double r = gsl_rng_uniform((gsl_rng *)rng); */

/*     // spin the wheel */
/*     double csum = 0.; */
/*     for (i=0; i<ivec->len; i++) { */
/* 	// get i^th index from ivec */
/* 	idx = g_array_index(ivec, int, i); */
/* 	// calculate probability for idx */
/* 	csum += vec[idx] / vec_sum; */
/* 	if (r < csum) { */
/* 	    // remove the selected index from the array of indices? */
/* 	    if (remove) */
/* 		g_array_remove_index_fast(ivec, i); */
/* 	    // return the index */
/* 	    return idx; */
/* 	} */
/*     } */

/*     // == none has been selected == */
/*     // remove last index from ivec? */
/*     if (remove) */
/* 	g_array_remove_index_fast(ivec, i-1); */
/*     // return last index */
/*     return idx; */
	
/* } */




/* // return a random int in [lo, hi) excluding the value EVAL */
/* int get_rnd_int_excluding_value(int lo, int hi, int eval, void *rng) { */
    
/*     int rnd = lo + gsl_rng_uniform_int((gsl_rng *)rng, hi - lo - 1); */

/*     return (rnd < eval ? rnd : rnd + 1); */
    
/* } */




/* // return a random int in [lo, hi) excluding the values in array EVALS */
/* // ** elements in EVALS should be in ascending order and mutually exclusive ** */
/* int get_rnd_int_excluding_values(int lo, int hi, int *evals, int evals_len, void *rng) { */

/*     int rnd = lo + gsl_rng_uniform_int((gsl_rng *)rng, hi - lo - evals_len); */
/*     int i; */

/*     for (i=0; i<evals_len; i++) { */
/* 	if (rnd < evals[i]) */
/* 	    break; */
/* 	rnd ++; */
/*     } */

/*     return rnd; */
    
/* } */



/* // return a double from a uniform distribution in [lo, hi] */
/* double get_rnd_double_with_bounds(double lo, double hi, void *rng) { */

/*     return lo + (hi-lo) * gsl_rng_uniform(rng); */
    
/* } */



