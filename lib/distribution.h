//==============================================
//       A DISTRIBUTION CLASS
//==============================================

typedef struct {

    GArray *_freq; // a vector of frequency recordings
    double _sval; // the starting value for computing the bins
    double _count; // total element count
    double dx; // the bin size

} dist_t;


// initialize a distribution
dist_t *dist_init(double sval, double dx);

// free a distribution
void dist_free(dist_t *dist);

// add a value to a distribution
void dist_add(dist_t *dist, double val);

// add a value to a distribution with a weight
void dist_add_w(dist_t *dist, double val, double weight);


// calculate the distribution's mean and variance
double dist_calc_mean(dist_t *dist, gboolean discrete);
double dist_calc_variance(dist_t *dist, double variance, gboolean discrete);

// save the distribution to a file with format:
//     line 0 : starting value
//     line 1 : dx
//     line 2 : tab-separated frequencies
// use Distribution.cload() from pylib/utils.py to parse file
void dist_dump(dist_t *dist, const char *fname);


// append the distribution as a string to STR
void dist_append_to_string(dist_t *dist, GString *str);


// empty the distribution 
void dist_reset(dist_t *dist);

