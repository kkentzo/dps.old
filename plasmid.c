#include <hdf5.h>

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <glib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include "../utils/clib.h"
#include "../utils/distribution.h"
#include "../utils/sampling.h"
#include "../utils/nvar.h"
#include "../utils/qhdf.h"

#include "pool.h"

#include "params.h"

#include "plasmid.h"

#include "cell.h"

#include "logger.h"









// ========================================================================
//                  A PLASMID REPLICATION PROFILE
//    characterized by a particular configuration of beta, alpha
// ========================================================================



// create a new plasmid replication profile
profile_t *profile_new(double beta, double kappa, double alpha, params_t *params)
{

    profile_t *profile = malloc(sizeof(profile_t));

    assert(profile);

    profile->beta = beta;
    profile->kappa = kappa;
    profile->alpha = alpha;
    
    profile->_ref_count = 1;

    profile->params = params;

    profile->parent = NULL;

    profile->step = params->step;

    // add profile to pool
    pool_add_profile(params->pool, profile);

    return profile;
    
}






// mutate the profile and return a new one
profile_t *profile_new_m(profile_t *profile) {

    int choice = g_array_index(profile->params->mwheel,
			       int,
			       gsl_rng_uniform_int(profile->params->rng,
						   profile->params->mwheel->len));

    double beta = profile->beta;
    double kappa = profile->kappa;
    double alpha = profile->alpha;

    switch (choice) 
	{

	case M_BETA:

	    beta = mutate(beta, profile->params->mut_rng, 0,
			  profile->params->max_beta, FALSE, profile->params->rng);
	    //((stats_t *)profile->params->stats)->mut_beta += 1;

	    break;

	case M_KAPPA:

	    kappa = mutate(kappa, profile->params->mut_rng, 0,
			   profile->params->max_kappa, FALSE, profile->params->rng);
	    //((stats_t *)profile->params->stats)->mut_beta += 1;

	    break;
	    
	case M_ALPHA:

	    alpha = mutate(alpha, profile->params->mut_rng, 0,
			   profile->params->max_alpha, FALSE, profile->params->rng);
	    //((stats_t *)profile->params->stats)->mut_alpha += 1;

	    break;
	    
	}


    // create mutant profile
    profile_t *mprofile = profile_new(beta, kappa, alpha, profile->params);
    // set the mutant's parent
    profile_set_parent(mprofile, profile);
    // and return mutant
    return mprofile;
	
}






// set the profile's parent (used in profile_new_m())
void profile_set_parent(profile_t *profile, profile_t *parent) {

    profile->parent = (struct profile_t *)parent;
    
}




// increase the object's reference count
profile_t *profile_retain(profile_t *profile) {

    // increase reference count
    profile->_ref_count ++;

    return profile;

}



// reduce the object's reference count
// if it goes to zero the memory is freed as well
profile_t *profile_release(profile_t *profile) {

    // reduce reference count
    if (-- profile->_ref_count == 0) {

	// remove entry from dictionary
	pool_remove_profile(profile->params->pool, profile);
	free(profile);
	profile = NULL;
	
    }

    return profile;
    
}






// increase/decrease the profile's copy number in the population
void profile_inc(profile_t *profile, int step) {

    pool_change_profile_count(profile->params->pool, profile, step);

}



void profile_dec(profile_t *profile, int step) {

    pool_change_profile_count(profile->params->pool, profile, - step);
    
}




// return the total cn of the profile
int profile_get_cn(profile_t *profile) {

    return pool_get_profile_cn(profile->params->pool, profile);
    
}





// append the profile's description in the supplied string
void profile_description(profile_t *profile, GString *st) {

    g_string_append_printf(st, "(%p %p %.3f %.3f %.3f %d %d %d)",
			   profile, (profile->parent ? profile->parent : 0),
			   profile->beta, profile->kappa, profile->alpha,
			   profile_get_cn(profile), profile->step,
			   profile->_ref_count);
    
}







// ========================================================================
//                      A PLASMID IN THE CELL
//             characterized by a profile and a copy number
// ========================================================================



// create a new plasmid with the supplied copy number
// ** retains the profile **
plasmid_t *plasmid_new(profile_t *profile, int cn) {

    plasmid_t *plasmid = malloc(sizeof(plasmid_t));

    // retain and store the profile
    plasmid->profile = profile_retain(profile);

    // set copy number
    plasmid->cn = cn;

    // increase profile's global copy number
    profile_inc(profile, cn);

    return plasmid;
    
}







// increase/decrease the plasmid's copy number
void plasmid_inc(plasmid_t *plasmid, int step) {

    plasmid->cn += step;
    profile_inc(plasmid->profile, step);

}



void plasmid_dec(plasmid_t *plasmid, int step) {

    plasmid->cn -= step;
    profile_dec(plasmid->profile, step);

}






// free the plasmid 
// ** release the profile **
void plasmid_free(plasmid_t *plasmid) {

    // decrease the profile's CN
    profile_dec(plasmid->profile, plasmid->cn);

    // release the profile
    profile_release(plasmid->profile);

    // free the plasmid
    free(plasmid);
    
}





// returns a GString with a description of the plasmid
// ** remember to free the string **
void plasmid_description(plasmid_t *plasmid, GString *st) {

    g_string_append_printf(st, "(%p %d)", plasmid->profile, plasmid->cn);
    
}



