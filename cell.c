#include <hdf5.h>

#include <string.h>
#include <assert.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>


#include "../utils/clib.h"
#include "../utils/sampling.h"
#include "../utils/distribution.h"
#include "../utils/bitvector.h"
#include "../utils/nvar.h"
#include "../utils/qhdf.h"


#include "pool.h"

#include "params.h"

#include "plasmid.h"

#include "cell.h"

#include "logger.h"

#include "population.h"









// =======================================================================
//                               CELL 
// =======================================================================


// create a new cell with the supplied parameters
cell_t *cell_new(GPtrArray *plasmids, params_t *params) {

    cell_t *cell = malloc(sizeof(cell_t));

    // store params
    cell->params = params;

    // store plasmids
    if (plasmids)
	cell->plasmids = plasmids;
    else
	cell->plasmids = g_ptr_array_new();

    // initialize array of mutants
    cell->mutants = g_ptr_array_new();

    // initialize cell descriptive statistics
    cell->stats = nvar_init(P_INTRA_IDX_ALL);

    // initialize other stuff
    cell_reset(cell);

    return cell;

}





// ====================================================================================
// reset the cell (to start a new cycle
void cell_reset(cell_t *cell) {


    // initialize internal state
    cell->omega = cell->params->omega;

    if (cell->params->omega_0_dev > 0)
	cell->omega_0 = cell->params->omega_0 +
	    gsl_ran_gaussian(cell->params->rng,
			     cell->params->omega_0_dev);
    else
	cell->omega_0 = cell->params->omega_0;

    cell->age = 0;

    cell->death = FALSE;
    cell->rdeath = FALSE;
    cell->division = FALSE;

    // DOMG will be calculated in the next cell_update()
    cell->domg = 0;

    // reset stats
    nvar_reset(cell->stats);

    // reset and recalculate copy number/sum_alpha
    cell->cn = 0;
    cell->sum_alpha = 0;
    plasmid_t *plasmid;
    int i;

    for (i=0; i<cell->plasmids->len; i++) {
    	// grab plasmid
    	plasmid = g_ptr_array_index(cell->plasmids, i);
	// update sums
	cell->cn += plasmid->cn;
	cell->sum_alpha += plasmid->cn * plasmid->profile->alpha;
	
    }

}




// ===============================================================================
// destroy the cell
void cell_free(cell_t *cell) {

    int i;    
    // free all plasmids 
    for (i=0; i<cell->plasmids->len; i++) 
	plasmid_free(g_ptr_array_index(cell->plasmids, i));

    // free the arrays themselves
    g_ptr_array_free(cell->plasmids, TRUE);
    g_ptr_array_free(cell->mutants, TRUE);

    // free the stats
    nvar_free(cell->stats);

    // free the cell
    free(cell);

}




// ===============================================================================
// add a plasmid in the cell of given profile with the supplied CN
// if the profile does not exist create a new plasmid
void cell_add_profile(cell_t *cell, profile_t *profile, int cn) {

    // update sums
    cell->cn += cn;
    cell->sum_alpha += cn * profile->alpha;

    // does the profile exist in the cell??
    int i;
    plasmid_t *plasmid;
    for (i=0; i<cell->plasmids->len; i++) {
	// grab plasmid
	plasmid = g_ptr_array_index(cell->plasmids, i);
	// is this the profile??
	if (profile == plasmid->profile) {
	    plasmid_inc(plasmid, cn);
	    return;
	}
    }

    // the profile does not exist -- add a new plasmid in both vectors!!
    g_ptr_array_add(cell->plasmids, plasmid_new(profile, cn));

}





// ===============================================================================
// select and return a plasmid randomly (prop. to n_i's)
plasmid_t *cell_select_plasmid(cell_t *cell) {

    if (cell->plasmids->len == 0)
	return NULL;

    plasmid_t *plasmid;
    double vals[cell->plasmids->len];
    int i;

    // determine probabilities
    for (i=0; i<cell->plasmids->len; i++) {
	// grab plasmid
	plasmid = g_ptr_array_index(cell->plasmids, i);
	// store copy number
	vals[i] = (double)plasmid->cn / cell->cn;

	assert(plasmid->cn > 0);
	
    }

    // select a random plasmid according to cals
    // store index in i
    sample(1, cell->plasmids->len, &i, vals, FALSE, FALSE, cell->params->rng);

    // return the plasmid
    return g_ptr_array_index(cell->plasmids, i);
    
}





// ===============================================================================
// divide the cell and return the daughter cell
cell_t *cell_divide(cell_t *cell) {

    int i, n_s;
    plasmid_t *plasmid;

    // create daughter cell
    cell_t *dcell = cell_new(NULL, cell->params);
        
    // whether to check for seg loss after division
    int parent_had_plasmids = cell->cn;

    // segregate all plasmids
    for (i=0; i<cell->plasmids->len; i++) {
	
    	// grab plasmid
    	plasmid = g_ptr_array_index(cell->plasmids, i);

	//assert(plasmid->cn > 0);

	// determine the copy number in the daughter cell
	if (cell->params->seg_type == SEG_BINOMIAL) 
	    n_s = gsl_ran_binomial(cell->params->rng, 0.5, plasmid->cn);
	else {
	    // split the plasmid in two 
	    n_s = plasmid->cn / 2;
	    // is there a remainder?? (odd number of plasmids)
	    if (plasmid->cn % 2 == 1)
		// decide randomly where the extra copy goes
		n_s += (gsl_rng_uniform(cell->params->rng) < 0.5 ? 0 : 1);
	}
	    
	// any plasmids leaving??
	if (n_s > 0)
	    // add segregated plasmids to daughter cell
	    cell_add_profile(dcell, plasmid->profile, n_s);
	
	// reduce copy number in the remaining plasmid
	plasmid_dec(plasmid, n_s);

	// are there any copies of the current plasmid
	// remaining in the cell??
	if (plasmid->cn == 0) {
	    
	    // remove plasmid from cell
	    g_ptr_array_remove_index_fast(cell->plasmids, i);
	    // update counter
	    -- i;
	    // free plasmid
	    plasmid_free(plasmid);
	    
	}
	
    }

    // reset current cell
    cell_reset(cell);

    // set the division flag (it will be reset in the next step's cell_grow())
    cell->division = TRUE;

    // check for segregation loss??
    if (parent_had_plasmids && (cell->cn == 0 || dcell->cn == 0))
	logger_register_events(cell->params->logger, E_LOSS, 1);

    // create and return daughter cell 
    return dcell;
    
}



// calculate cellular growth and set the DIV/DEATH flags
// if the cell divides returns a pointer
// to the newly created daughter cell, otherwise CELL is returned
cell_t *cell_grow(cell_t *cell) {

    // calculate new delta omega
    cell->domg = cell->omega_0 +			\
	cell->params->phi * cell->cn / (cell->params->lambda + cell->cn) - \
	cell->params->gamma * cell->cn - cell->params->gamma_alpha * cell->sum_alpha;

    // update omega
    cell->omega += cell->domg;

    // === update flags ===
    // do we have a division/death event??
    cell->division = (cell->omega >= 2 * cell->params->omega);

    cell->death = (cell->omega < 0 || cell->cn > cell->params->max_cn);

    cell->rdeath = FALSE;

    // === register division or death ===
    if (cell->death) 
	logger_register_events(cell->params->logger, E_DEATH, 1);
    else if (cell->division) 
	// register cell division
	logger_register_division(cell->params->logger, cell);

    if (cell->division)
	assert(! cell->death);

    // perform division???
    return (cell->division ? cell_divide(cell) : cell);

}






// update the cell's state (plasmid replication/conjugation)
// add conjugant plasmids into conjugation pool
// register cell in global/intra/inter stats
void cell_update(cell_t *cell, pool_t *cpool, gboolean rdeath) {

    // update the RDEATH flag
    cell->rdeath = rdeath;

    // perform conjugation (only if the cell is not dead or dividing)
    plasmid_t *conjugant = NULL;

    if ((! cell->death) && cell->cn > 0 &&
	cell->params->pconj > 0 &&
	gsl_rng_uniform(cell->params->rng) < cell->params->pconj)
	{

	    // select a random plasmid from the cell
	    conjugant = cell_select_plasmid(cell);
	    // add its profile to the conjugation pool
	    pool_change_profile_count(cpool, conjugant->profile, 1);
	    
	} 

    
    // store sum_alpha locally - we want the same sum_alpha
    // for the replication of all profiles
    double sum_alpha = cell->sum_alpha;
    int cn;
    // store the total extra CN (for updating cell->cn at the end)
    // **ATTENTION** cell->cn does (and should) not change during plasmid replication
    cell->total_extra_cn = 0;

    plasmid_t *plasmid;
    int i, j;
    double r;
    double n_r, n_m;
    profile_t *profile;
    profile_t *mut_profile;

    // reset intra-cellular stats (they will be recalculated)
    nvar_reset(cell->stats);

    double pvals[P_INTRA_IDX_ALL];
    
    // consider all plasmids
    for (i=0; i<cell->plasmids->len; i++) {
	
	// grab plasmid
	plasmid = g_ptr_array_index(cell->plasmids, i);

	// store plasmid cn
	cn = plasmid->cn;

	// calculate per-plasmid replication rate
	r = plasmid->profile->beta /				\
	    (1. + plasmid->profile->kappa * sum_alpha /
	     (cell->params->dilution ? cell->omega : 1.));

	// calculate number of replication events for the plasmid type (i.e. R times CN)
	n_r = gsl_ran_poisson(cell->params->rng, cn * r);

	// register profile in intra-cellular stats
	pvals[P_INTRA_IDX_CN] = cn;
	pvals[P_INTRA_IDX_BETA] = plasmid->profile->beta;
	pvals[P_INTRA_IDX_KAPPA] = plasmid->profile->kappa;
	pvals[P_INTRA_IDX_ALPHA] = plasmid->profile->alpha;
	pvals[P_INTRA_IDX_RATE] = r;
	pvals[P_INTRA_IDX_BK] = plasmid->profile->beta / plasmid->profile->kappa;
	pvals[P_INTRA_IDX_N_R] = (double)n_r / cn;
	pvals[P_INTRA_IDX_HT] = (plasmid == conjugant ? 1. / cn : 0);
	pvals[P_INTRA_IDX_DEATH] = (cell->death || cell->rdeath ?
				    1. + pvals[P_INTRA_IDX_N_R] :
				    0);

	// this represents the per-plasmid profile's
	// representation at the next time step
	pvals[P_INTRA_IDX_FITNESS] =
	    + 1.
	    + pvals[P_INTRA_IDX_N_R]
	    + pvals[P_INTRA_IDX_HT]
	    - pvals[P_INTRA_IDX_DEATH];

	// transmission biases are of interest only at the GLOBAL level
	// do not bother with them here
	pvals[P_INTRA_IDX_TBETA] = 0;
	pvals[P_INTRA_IDX_TKAPPA] = 0;
	pvals[P_INTRA_IDX_TALPHA] = 0;
	pvals[P_INTRA_IDX_TBK] = 0;

	// update total number of replications
	cell->total_extra_cn += n_r;

	// update plasmid's copy number
	plasmid_inc(plasmid, n_r);

	// register replications
	logger_register_events(cell->params->logger, E_REP, n_r);

	// reset mutation counter
	n_m = 0;

	// do we consider mutations (bka)?? (only if the cell does not die)
	if (!(cell->death || cell->rdeath) && n_r > 0 && cell->params->mwheel->len > 0) {
    
	    // out of all the replication events
	    // how many are the mutation events??
	    n_m = gsl_ran_binomial(cell->params->rng, cell->params->mu, n_r);

	    // decrease the copy number of the plasmid by n_m
	    plasmid_dec(plasmid, n_m);
	    // and decrease the extra cn counter (because the mutant is
	    // going to be added to cell->cn in cell_add_profile)
	    cell->total_extra_cn -= n_m;

	    // register mutation events
	    if (n_m > 0)
		logger_register_events(cell->params->logger, E_MUT, n_r);

	    // create n_m mutants of the current profile and add them to the mutants array
	    // remember to release profile when emptying the mutant array
	    for (j=0; j<n_m; j++) {
		mut_profile = profile_new_m(plasmid->profile);
		// register mutation globally
		pool_register_mutation(cell->params->pool,
				       plasmid->profile,
				       mut_profile);
		// add the mutant to array of mutants 
		g_ptr_array_add(cell->mutants, mut_profile);
	    }

	}

	// register profile in the pool
	pool_register_profile(cell->params->pool, plasmid->profile, cn,
			      pvals[P_INTRA_IDX_FITNESS], n_m);

	// add profile to INTRA stats
	nvar_add(cell->stats, pvals, cn);

	// update global sum_alpha and sum_mobility (these will be updated separately for
	// potential mutants in cell_add_profile())
	cell->sum_alpha += (n_r - n_m) * plasmid->profile->alpha;

    }

    // register the intra-cellular stats in the pool
    pool_register_cell_stats(cell->params->pool, cell->stats);

    // register cell in stats
    logger_register_cell(cell->params->logger, cell);

    // update copy number (now that the cell is registered in logger)
    cell->cn += cell->total_extra_cn;

    // add mutant profiles to cell
    for (i=0; i<cell->mutants->len; i++) {
	// grab mutant profile
	profile = g_ptr_array_index(cell->mutants, i);
	// add profile to host
	cell_add_profile(cell, profile, 1);
	// release profile (it's now owned by a plasmid)
	profile_release(profile);
    }

    // empty array of mutants
    g_ptr_array_set_size(cell->mutants, 0);

    // update cell's age
    cell->age ++ ;

}





// append a description of the cell to gstring
// form: (p1_descr p2_descr ...)
void cell_description(cell_t *cell, GString *st) {

    int i;
    // write cell params
    g_string_append_printf(st, "(%.3e %.3e %d %d %d ", cell->omega, cell->omega_0,
			   cell->age, (cell->division ? 1 : 0),
			   (cell->death ? 1 : 0));
    //g_string_append_c(st, '(');
    g_string_append_printf(st, " (");

    // write plasmids
    for (i=0; i<cell->plasmids->len; i++) {
	plasmid_description(g_ptr_array_index(cell->plasmids, i), st);
	if (i < cell->plasmids->len - 1)
	    g_string_append_c(st, ' ');	
    }
    
    g_string_append_printf(st, "))");
    
}
