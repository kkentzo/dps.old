#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <hdf5.h>
#include <glib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include "../utils/distribution.h"
#include "../utils/nvar.h"
#include "../utils/qhdf.h"
#include "../utils/clib.h"


#include "pool.h"

#include "params.h"

#include "plasmid.h"

#include "cell.h"

#include "logger.h"




// initialize the pool
pool_t *pool_new() {

    pool_t *pool = malloc(sizeof(pool_t));

    // initialize HASH TABLE
    pool->profiles = g_hash_table_new_full(g_direct_hash,
					   g_direct_equal,
					   NULL,
					   free);

    pool->pid = 0;

    pool->cn = 0;

    // initialize matrices
    pool->stats = nvar_init(P_INTRA_IDX_ALL);

    pool->size = 0;

    // initialize buffers
    pool->M = calloc(P_INTRA_IDX_ALL, sizeof(double));
    pool->V = calloc(P_INTRA_IDX_ALL, sizeof(double));
    pool->C = calloc(P_INTRA_IDX_ALL*(P_INTRA_IDX_ALL-1)/2, sizeof(double));

    return pool;
    
}



// free the pool
void pool_free(pool_t *pool) {

    // free the hash table
    g_hash_table_destroy(pool->profiles);

    // free the stats
    nvar_free(pool->stats);

    // free the buffers
    free(pool->M);
    free(pool->V);
    free(pool->C);

    // free the struct
    free(pool);
    
}



// empty the pool
void pool_empty(pool_t *pool) {

    g_hash_table_remove_all(pool->profiles);

    nvar_reset(pool->stats);

    pool->cn = 0;    

    pool->size = 0;

}




// add a profile to the pool (with CN=0)
void pool_add_profile(pool_t *pool, void *profile) {

    pdata_t *pdata = g_hash_table_lookup(pool->profiles, profile);
    
    if (pdata == NULL) {
	pdata = malloc(sizeof(pdata_t));

	// set profile pid and increase global pid
	pdata->pid = pool->pid ++;
	pdata->cn = 0;

	pdata->mbeta_sum = 0;
	pdata->mkappa_sum = 0;
	pdata->malpha_sum = 0;
	pdata->mbk_sum = 0;
	pdata->fitness_sum = 0;

	pdata->pcn = 0;
	pdata->mutants = 0;

	/* pdata->include = FALSE; */

	g_hash_table_insert(pool->profiles, profile, pdata);

	pool->size += 1;
    }
    
}



// remove a profile 
void pool_remove_profile(pool_t *pool, void *profile) {

    int cn = pool_get_profile_cn(pool, profile);

    pool_change_profile_count(pool, profile, - cn);

    g_hash_table_remove(pool->profiles, profile);

    pool->cn -= cn;

    pool->size -= 1;
    
}




// increase/decrease the copy number of a profile
// the profile is added to the pool if it's not there already
void pool_change_profile_count(pool_t *pool, void *profile, int step) {

    pdata_t *pdata = g_hash_table_lookup(pool->profiles, profile);
    // add the profile to the pool if necessary
    if (pdata == NULL) {
	pool_add_profile(pool, profile);
	pdata = g_hash_table_lookup(pool->profiles, profile);
    }

    pdata->cn += step;

    pool->cn += step;

}





// return the CN of a profile
int pool_get_profile_cn(pool_t *pool, void *profile) {

    pdata_t *pdata = g_hash_table_lookup(pool->profiles, profile);

    return (pdata ? pdata->cn : 0);
    
}





// register intra-cellular stats
void pool_register_cell_stats(pool_t *pool, nvar_t *stats) {

    nvar_add_stats(pool->stats, stats);
    
}




// register replication (for computing transmission biases)
// ** offspring should include any conjugation events **
void pool_register_profile(pool_t *pool, void *profile, int cn,
			   double fitness, int mutants)
{

    // grab profile's pdata
    pdata_t *pdata = g_hash_table_lookup(pool->profiles, profile);
    assert(pdata);

    pdata->fitness_sum += cn * fitness;
    pdata->pcn += cn;
    pdata->mutants += mutants;
    
}



// register mutation (for computing transmission biases)
void pool_register_mutation(pool_t *pool, void *wt_profile, void *mut_profile) {

    // grab profile's pdata
    pdata_t *pdata = g_hash_table_lookup(pool->profiles, wt_profile);

    assert(pdata);

    pdata->mbeta_sum += ((profile_t *)mut_profile)->beta;
    pdata->mkappa_sum += ((profile_t *)mut_profile)->kappa;
    pdata->malpha_sum += ((profile_t *)mut_profile)->alpha;
    pdata->mbk_sum += (((profile_t *)mut_profile)->beta /
		       ((profile_t *)mut_profile)->kappa);

}





void pool_description(pool_t *pool, GString *st) {

    GHashTableIter iter;
    pdata_t *pdata;
    profile_t *profile;

    g_hash_table_iter_init(&iter, pool->profiles);

    g_string_assign(st, "(");

    while (g_hash_table_iter_next(&iter, (void **)&profile, (void **)&pdata)) {

	assert(profile);

	profile_description(profile, st);
	g_string_append_c(st, ' ');
	
    }

    g_string_append_c(st, ')');
    
}



// reset ht's and vt's for all profiles
void pool_update(pool_t *pool, int step, void *_params) {


    params_t *params = _params;
    GHashTableIter iter;
    profile_t *profile;
    pdata_t *pdata;

    g_hash_table_iter_init(&iter, pool->profiles);

    int freqs[params->contenders];
    memset(freqs, 0, params->contenders * sizeof(int));

    // initialize stats for the transmission biases
    nvar_t *tbiases = nvar_init(4);
    double pvals[4];
    double f;

    while (g_hash_table_iter_next(&iter, (void **)&profile, (void **)&pdata)) {

    	// == update joint distributions in logger ==
	hdf_histogram_add2(&((logger_t *)params->logger)->hist_beta_kappa,
			   profile->beta, profile->kappa, pdata->cn);
	hdf_histogram_add2(&((logger_t *)params->logger)->hist_beta_alpha,
			   profile->beta, profile->alpha, pdata->cn);
	hdf_histogram_add2(&((logger_t *)params->logger)->hist_kappa_alpha,
			   profile->kappa, profile->alpha, pdata->cn);

	// update frequencies (in competition mode)
	if (params->compete) {

	    assert(pdata->pid < params->contenders);
	    freqs[pdata->pid] = pdata->cn;
	    
	}

	// update transmission biases (weight by number of offspring)
	if (pdata->fitness_sum > 0) {
	    // calculate average profile fitness
	    f = pdata->fitness_sum / pdata->pcn;
	    // calculate the profile's contribution to the overall sum of the tbias

	    pvals[0] = ((pdata->mbeta_sum +
			 (pdata->fitness_sum - pdata->mutants) * profile->beta) / \
			pdata->fitness_sum - profile->beta) * f;

	    pvals[1] = ((pdata->mkappa_sum +
			 (pdata->fitness_sum - pdata->mutants) * profile->kappa) / \
			pdata->fitness_sum - profile->kappa) * f;

	    pvals[2] = ((pdata->malpha_sum +
			 (pdata->fitness_sum - pdata->mutants) * profile->alpha) / \
			pdata->fitness_sum - profile->alpha) * f;

	    pvals[3] = ((pdata->mbk_sum +
			 (pdata->fitness_sum - pdata->mutants) *	\
			 profile->beta / profile->kappa) /		\
			pdata->fitness_sum - profile->beta / profile->kappa) * f;

	    // add pvals to tbiases (weighted by number of offspring)
	    nvar_add(tbiases, pvals, pdata->fitness_sum);

	} 

	
	/* // reset profile's pdata */
	pdata->mbeta_sum = 0;
	pdata->mkappa_sum = 0;
	pdata->malpha_sum = 0;
	pdata->mbk_sum = 0;
	pdata->fitness_sum = 0;

	pdata->pcn = 0;
	pdata->mutants = 0;

	/* pdata->include = FALSE; */

    }

    // set the mean beta, kappa and alpha in logger (for screen output only)
    ((logger_t *)params->logger)->beta = nvar_calc_mean(pool->stats, P_INTRA_IDX_BETA);
    ((logger_t *)params->logger)->kappa = nvar_calc_mean(pool->stats, P_INTRA_IDX_KAPPA);
    ((logger_t *)params->logger)->alpha = nvar_calc_mean(pool->stats, P_INTRA_IDX_ALPHA);

    // write the stats??
    if (step % params->log_every == 0) {
	//printf("G=%.10e | ", nvar_calc_covariance(pool->stats, P_INTRA_IDX_BETA, P_INTRA_IDX_FITNESS, FALSE));
	// write the stats to buffers
	nvar_get_statistics(pool->stats, pool->M, pool->V, pool->C, FALSE, 1);
	// replace transmission biases
	pool->M[P_INTRA_IDX_TBETA] = nvar_calc_mean(tbiases, 0);
	pool->M[P_INTRA_IDX_TKAPPA] = nvar_calc_mean(tbiases, 1);
	pool->M[P_INTRA_IDX_TALPHA] = nvar_calc_mean(tbiases, 2);
	pool->M[P_INTRA_IDX_TBK] = nvar_calc_mean(tbiases, 3);
	
	// write the buffers to the tables
	hdf_table_append_record(&((logger_t *)params->logger)->tbl_global_m, pool->M);
	hdf_table_append_record(&((logger_t *)params->logger)->tbl_global_v, pool->V);
	hdf_table_append_record(&((logger_t *)params->logger)->tbl_global_c, pool->C);
	// update the competition table in logger
	if (params->compete) 
	    hdf_table_append_record(&((logger_t *)params->logger)->tbl_competition, freqs);
	// reset the stats
	nvar_reset(pool->stats);

    }

    // free the tbiases
    nvar_free(tbiases);

}








// distribute the conjugants among recipients
void pool_distribute(pool_t *pool, GPtrArray *recipients, void *_params) {

    // cast and store params
    params_t *params = _params;

    GHashTableIter iter;
    profile_t *profile;
    pdata_t *pdata;
    int i, idx;
    cell_t *recipient;
    int events = 0;

    g_hash_table_iter_init(&iter, pool->profiles);

    // for each profile in the conjugation pool
    while (g_hash_table_iter_next(&iter, (void **)&profile, (void **)&pdata)) {

	assert(pdata->cn > 0);

	// perform transfers
	for (i=0; i<pdata->cn; i++) {

	    // choose a random recipient cell
	    idx = gsl_rng_uniform_int(params->rng, recipients->len);
	    recipient = g_ptr_array_index(recipients, idx);
	    // transmit the profile
	    cell_add_profile(recipient, profile, 1);
	    
	}

	// update event counter
	events += pdata->cn;

    }

    logger_register_events(params->logger, E_HT, events);
    
}
