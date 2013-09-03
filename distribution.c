#include <stdlib.h>
#include <stdio.h>
#include <glib.h>


#include "distribution.h"


// initialize a distribution
dist_t *dist_init(double sval, double dx) {

    // claim mem for object
    dist_t *dist = malloc(sizeof(dist_t));

    // initialize vector of frequencies
    dist->_freq = g_array_new(FALSE, FALSE, sizeof(double));
    // initialize starting value
    dist->_sval = sval;
    // store dx
    dist->dx = dx;
    // initialize total count
    dist->_count = 0;

    return dist;
    
}


// free a distribution
void dist_free(dist_t *dist) {

    if (dist) {
	// free the vector
	g_array_free(dist->_freq, TRUE);
	// free the dist
	free(dist);
    }
    
}




// val is before the start of the distribution
void _dist_prepend_(dist_t *dist, double val, double weight) {

    double zero = 0;

    while (1) {
	// create a new bin on the left
	g_array_prepend_val(dist->_freq, zero);
	// adjust sval
	dist->_sval -= dist->dx;
	// do we have the correct bin??
	if (val >= dist->_sval)
	    break; // got it!
    }	

    // record the element by increasing the frequency
    // of the last added bin
    g_array_index(dist->_freq, double, 0) += weight;

}




// val is after the end of the distribution
void _dist_append_(dist_t *dist, double val, double weight) {

    double zero = 0;

    // keep appending bins until val falls inside one of them
    while (1) {
	// create a new bin on the right
	g_array_append_val(dist->_freq, zero);
	// does val fall inside??
	if (val < dist->_sval + dist->_freq->len * dist->dx)
	    break; // got it!
    }	

    // record the element by increasing the frequency
    // of the last added bin
    g_array_index(dist->_freq, double, dist->_freq->len - 1) += weight;
    
}




// add a value to a distribution
void dist_add(dist_t *dist, double val) {

    dist_add_w(dist, val, 1);

}





// add a value to a distribution
void dist_add_w(dist_t *dist, double val, double weight) {

    int i;

    if (val < dist->_sval)
	_dist_prepend_(dist, val, weight);
    else if (val > dist->_sval + dist->_freq->len * dist->dx)
	_dist_append_(dist, val, weight);
    else 
	// add the value in the distribution
	for (i=0; i<dist->_freq->len; i++) {
	    if (val < dist->_sval + (i+1) * dist->dx) {
		// got the bin : increase frequency of element by one
		g_array_index(dist->_freq, double, i) += weight;
		break;
	    }
	}


    // increase the total sum of elements
    dist->_count += weight;

}




// calculate the distribution's mean
double dist_calc_mean(dist_t *dist, gboolean discrete) {

    int i;

    double mean = 0;
    double val;

    for (i=0; i<dist->_freq->len; i++) {
	if (discrete)
	    val = (dist->_sval + i*dist->dx);
	else
	    val = (dist->_sval + i*dist->dx + dist->dx / 2.);
	mean += val * g_array_index(dist->_freq, double, i) / dist->_count;
    }

    return mean;
}




// calculate the distribution's variance (supply the mean)
double dist_calc_variance(dist_t *dist, double mean, gboolean discrete) {

    if (dist->_count < 1e-10)
	return 0;

    int i;

    double variance = 0;
    double val;
    double prob;

    for (i=0; i<dist->_freq->len; i++) {
	if (discrete)
	    val = (dist->_sval + i*dist->dx);
	else
	    val = (dist->_sval + i*dist->dx + dist->dx / 2.);
	prob = g_array_index(dist->_freq, double, i) / dist->_count;
	variance += (val - mean) * (val - mean) * prob;
	
    }

    return variance;
    
}


		 

// save the distribution to a file
void dist_dump(dist_t *dist, const char *fname) {

    FILE *f = fopen(fname, "w");

    // write sval
    fprintf(f, "%f\n", dist->_sval);
    // write dx
    fprintf(f, "%f\n", dist->dx);
    // write frequencies
    int i;
    for (i=0; i<dist->_freq->len; i++)
	fprintf(f, "%.1f\t", g_array_index(dist->_freq, double, i));
    
    fclose(f);
    
}



// append the distribution as a string to STR
void dist_append_to_string(dist_t *dist, GString *str) {

    // append sval and dx
    g_string_append_printf(str, "%f:%f:", dist->_sval, dist->dx);
    // append frequencies
    int i;
    for (i=0; i<dist->_freq->len; i++)
	g_string_append_printf(str, "%.1f\t", g_array_index(dist->_freq, double, i));
    
}




// empty the distribution 
void dist_reset(dist_t *dist) {

    g_array_set_size(dist->_freq, 0);
    dist->_count = 0;
    
}
