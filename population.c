#include <hdf5.h>
#include <hdf5_hl.h>

#include <time.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <glib.h>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include "../utils/sampling.h"
#include "../utils/distribution.h"
#include "../utils/clib.h"
#include "../utils/bitvector.h"
#include "../utils/nvar.h"
#include "../utils/qhdf.h"

#include "pool.h"

#include "params.h"

#include "plasmid.h"
#include "cell.h"

#include "logger.h"

#include "population.h"







void run(params_t *params) {

    cell_t *cell, *dcell;

    int i, rdeath;

    pool_t *cpool = pool_new();

    // start counting the time
    params->start = time(NULL);

    int the_dead[2 * params->psize]; // indices of random-death hosts
    bitvector_t *dvec = bv_new(2 * params->psize);

    GPtrArray *dead_cells = g_ptr_array_new();

    // run steps
    for (params->step=0; params->step<params->steps; params->step++) {

	if (global_abort)
	    break;

	// update omegas and count number of division and death events
	// perform cell divisions and cell deaths (**not** random deaths)
	for (i=0, params->ndiv=0, params->ndeath=0; i<params->cells->len; i++)
	    {
		// grab cell
		cell = g_ptr_array_index(params->cells, i);
		// grow the cell
		dcell = cell_grow(cell);
		// do we have a division event??
		if (dcell != cell) {
		    g_ptr_array_add(params->dcells, dcell);
		    ++ params->ndiv;
		}

		// do we have a death event??
		if (cell->death) {
		    ++ params->ndeath;
		    g_ptr_array_add(dead_cells,
				    g_ptr_array_remove_index_fast(params->cells, i));
		    -- i;
		}
	    }

	// append all daughter cells to params->cells
	for (i=0; i<params->dcells->len; i++)
	    g_ptr_array_add(params->cells, g_ptr_array_index(params->dcells, i));
	g_ptr_array_set_size(params->dcells, 0);


	// ** At this point params->cells does not contain dead cells **

	// calculate number of random cell deaths
	// (i.e. determine number of excess cells)
	rdeath = (params->cells->len > params->psize ?
		  params->cells->len % params->psize :
		  0);

	// select cells to be randomly killed (includes daughter cells)
	sample(rdeath, params->cells->len, the_dead, NULL, FALSE, FALSE, params->rng);
	// reset the bitvector
	bv_set_all_off(dvec);
	// set the bitvector bits on for the cells to be randomly killed
	for (i=0; i<rdeath; i++) 
	    bv_set_on(dvec, the_dead[i]);

	// === REGISTER THE GLOBAL/INTER/INTRA STATS FOR ALL HOSTS ===
	
	// === update all dead cells  ===
	for (i=0; i<dead_cells->len; i++) {
	    // grab cell
	    cell = g_ptr_array_index(dead_cells, i);
	    // update cell
	    cell_update(cell, cpool, FALSE);
	}

	// === update all living cells ===
	for (i=0; i<params->cells->len; i++) {
	    // grab cell
	    cell = g_ptr_array_index(params->cells, i);
	    // update cell
	    cell_update(cell, cpool, bv_is_on(dvec, i));
	}

	// do we have conjugant plasmids??
	if (cpool->size > 0) {
	    // distribute conjugants
	    pool_distribute(cpool, params->cells, params);
	    // empty cpool
	    pool_empty(cpool);
	}

	// free the cells that have randomly died
	for (i=0; i<params->cells->len; i++) {
	    // grab cell
	    cell = g_ptr_array_index(params->cells, i);
	    if (cell->rdeath) {
		cell_free(cell);
		g_ptr_array_remove_index_fast(params->cells, i);
		--i;
	    }
	}

	// free the dead cells
	for (i=0; i<dead_cells->len; i++)
	    cell_free(g_ptr_array_index(dead_cells, i));
	g_ptr_array_set_size(dead_cells, 0);

	// is population extinct??
	if (params->cells->len == 0) {
	    printf("step %d | Population is exinct!\n", params->step);
	    break;
	}

	// are plasmid extinct??
	if ((int)params->pool->cn == 0) {
	    printf("step %d | Plasmids are exinct!\n", params->step);
	    break;
	}

	// is competition over??
	if (params->compete && params->pool->size == 1) {
	    printf("step %d | Competition is over\n", params->step);
	    global_abort = TRUE;
	}

	// update the histograms and log the global statistics
	pool_update(params->pool, params->step, params);

	// log statistics
	logger_log(params->logger, params->step);

    }

    // free the bitvector
    bv_free(dvec);

    // free the conjugants' pool
    pool_free(cpool);

    // free the dead_cells array
    g_ptr_array_free(dead_cells, TRUE);

    // calculate duration
    params->duration = difftime(time(NULL), params->start);

    // save the population
    save_population(params, NULL);

}















// ==================================================
//       Population Initialization and Saving
// ==================================================




// assemble the profile and hosts strings and call logger_save_population()
void save_population(params_t *params, const char *fname) {

    // population string
    GString *pop_st = g_string_new("(");
    // profiles string
    //GString *prof_st = g_string_new("(");
    //GHashTable *profiles = g_hash_table_new(g_direct_hash, g_direct_equal);
    //profile_t *profile;
    cell_t *cell;

    int i; //j;

    // for all cells within the thread
    for (i=0; i<params->cells->len; i++) {
	// grab cell
	cell = g_ptr_array_index(params->cells, i);
	// get cell's description
	cell_description(cell, pop_st);
	g_string_append_c(pop_st, ' ');
	// for all plasmids within the cell
	/* for (j=0; j<cell->plasmids->len; j++) { */
	/*     // get plasmid profile */
	/*     profile = ((plasmid_t *)g_ptr_array_index(cell->plasmids, j))->profile; */
	/*     if (g_hash_table_lookup(profiles, profile) == NULL) */
	/* 	g_hash_table_insert(profiles, profile, NULL); */
	/* } */
    }

    g_string_append_c(pop_st, ')');

    /* GHashTableIter iter; */
    /* g_hash_table_iter_init(&iter, profiles); */
    /* void *val; */
    /* while (g_hash_table_iter_next(&iter, (void **)&profile, &val)) { */
    /* 	profile_description(profile, prof_st); */
    /* 	g_string_append_c(prof_st, ' '); */
    /* } */

    /* g_string_append_c(prof_st, ')'); */

    /* // free the hash table */
    /* g_hash_table_destroy(profiles); */


    GString *prof_st = g_string_new("");

    pool_description(params->pool, prof_st);

    // where to save??
    if (fname) {

	// write to file
	FILE *fp = fopen(fname, "w");
	fprintf(fp, "%s\n%s", prof_st->str, pop_st->str);
	fclose(fp);
	
    } else {
	
	// pass the strings to logger
	logger_save_population(params->logger, prof_st, pop_st);
	
    }

    // free the strings
    g_string_free(prof_st, TRUE);
    g_string_free(pop_st, TRUE);
    
}









// ========================================================================
// load the population from params->load_from
// and add all cells in CELLS
void load_population_from_file(params_t *params) {

    // open HDF file
    hid_t file_id = H5Fopen(params->load_from->str, H5F_ACC_RDONLY, H5P_DEFAULT);

    hid_t dataset_id, filetype_id;
    size_t sdim;
    char *buffer;

    GString *profiles_string, *hosts_string;


    // ==== READ IN PROFILE STRING =====
    // open profile dataset
    dataset_id = H5Dopen1(file_id, "/population/profiles");
    // get the datatype and its size
    filetype_id = H5Dget_type(dataset_id);
    sdim = H5Tget_size(filetype_id);
    // close filetype and dataset
    H5Tclose(filetype_id);
    H5Dclose(dataset_id);
    // initialize profile string
    buffer = malloc(sdim * sizeof(char));
    // read in string
    H5LTread_dataset_string(file_id, "/population/profiles", buffer);
    // create GString for profiles
    profiles_string = g_string_new(buffer);
    // free buffer
    free(buffer);


    // ==== READ IN HOSTS STRING =====
    // open profile dataset
    dataset_id = H5Dopen1(file_id, "/population/hosts");
    // get the datatype and its size
    filetype_id = H5Dget_type(dataset_id);
    sdim = H5Tget_size(filetype_id);
    // close filetype and dataset
    H5Tclose(filetype_id);
    H5Dclose(dataset_id);
    // initialize profile string
    buffer = malloc(sdim * sizeof(char));
    // read in string
    H5LTread_dataset_string(file_id, "/population/hosts", buffer);
    // create GString for profiles
    hosts_string = g_string_new(buffer);
    // free buffer
    free(buffer);

    //printf("profiles=%s\nhosts=%s\n", profiles_string->str, hosts_string->str);

    // close file
    H5Fclose(file_id);

    // read in prof_st and construct profiles
    GHashTable *profiles = g_hash_table_new_full((GHashFunc)g_string_hash,
						 (GEqualFunc)g_string_equal,
						 (GDestroyNotify)g_string_free_func,
						 (GDestroyNotify)profile_release);

    int i, j;
    GPtrArray *tokens = g_string_tokenize(profiles_string);
    if (tokens) {

	GPtrArray *ptokens, *qtokens, *vtokens;
	GString *key;
	double beta, kappa, alpha;
	for (i=0; i<tokens->len; i++) {
	    // tokenize string
	    ptokens = g_string_tokenize(g_ptr_array_index(tokens, i));
	    key = g_ptr_array_index(ptokens, 0);

	    // ATTENTION : token #1 is the pointer to the profile's parent
	    //             which does not concern us right now
	    
	    if (params->f_beta)
		beta = params->beta;
	    else
		beta = atof(((GString *)g_ptr_array_index(ptokens, 2))->str);

	    if (params->f_kappa)
		kappa = params->kappa;
	    else
		kappa = atof(((GString *)g_ptr_array_index(ptokens, 3))->str);
	    
	    if (params->f_alpha)
		alpha = params->alpha;
	    else
		alpha = atof(((GString *)g_ptr_array_index(ptokens, 4))->str);

	    
	    if (tokens->len < 5) // why 5??? ==> aah, do not print more than 5 profiles, OK.
		printf("Adding profile : beta=%.3f | kappa=%.3f | alpha=%.3f \n",
		       beta, kappa, alpha);

	    // construct profile and add it to hash table
	    g_hash_table_insert(profiles,
				g_string_new(key->str),
				profile_new(beta, kappa, alpha, params));

	    // free ptokens
	    g_ptr_array_free_with_func(ptokens, g_string_free_func);

	}

	g_ptr_array_free_with_func(tokens, g_string_free_func);

	
	// === now load the population ===
	cell_t *cell;
	profile_t *profile;
	int cn;
	GString *st;
    
	tokens = g_string_tokenize(hosts_string);
	// this shouldn't be NULL
	assert(tokens);

	// for each cell
	for (i=0; i<tokens->len; i++) {

	    // tokenize cell contents
	    ptokens = g_string_tokenize(g_ptr_array_index(tokens, i));
	    assert(ptokens);

	    // remove empty tokens (this was a bug that prevented the
	    // parsing of plasmids into cells)
	    for (j=0; j<ptokens->len; j++) {
		// get j^th token
		st = g_ptr_array_index(ptokens, j);
		// is string empty??
		if (st->len == 0) {
		    g_string_free(st, TRUE);
		    g_ptr_array_remove_index(ptokens, j);
		    --j;
		} /* else { */
		/*     printf("|%d : _%s_ ", j, st->str); */
		/*     printf("\n"); */
		/* } */
	    }

	    // create cell
	    cell = cell_new(NULL, params);

	    // get value of omega and omega_0 and save to cell
	    cell->omega = atof(((GString *)g_ptr_array_index(ptokens, 0))->str);
	    cell->omega_0 = atof(((GString *)g_ptr_array_index(ptokens, 1))->str);
	    cell->age = atoi(((GString *)g_ptr_array_index(ptokens, 2))->str);

	    // !! skip division and death tokens (3 and 4)

	    // tokenize plasmids
	    vtokens = g_string_tokenize(g_ptr_array_index(ptokens, 5));

	    if (vtokens == NULL) {
		// it's plasmid free
		g_ptr_array_add(params->cells, cell);
		continue;
	    }

	    // for each plasmid in the cell
	    for (j=0; j<vtokens->len; j++) {

		qtokens = g_string_tokenize(g_ptr_array_index(vtokens, j));
		// get profile key
		key = g_ptr_array_index(qtokens, 0);
		cn = atoi(((GString *)g_ptr_array_index(qtokens, 1))->str);
		// retrieve profile from hash table
		profile = g_hash_table_lookup(profiles, key);
		assert(profile);
		// create plasmid and add it to cell
		cell_add_profile(cell, profile, cn);
		// free qtokens
		g_ptr_array_free_with_func(qtokens, g_string_free_func);
		
	    }

	    // add cell to cells
	    g_ptr_array_add(params->cells, cell);

	    // free vtokens
	    g_ptr_array_free_with_func(vtokens, g_string_free_func);
	    // free ptokens
	    g_ptr_array_free_with_func(ptokens, g_string_free_func);
	
	}

	g_ptr_array_free_with_func(tokens, g_string_free_func);	

    }


    // free hash table
    g_hash_table_destroy(profiles);

    // free stuff
    g_string_free(hosts_string, TRUE);
    g_string_free(profiles_string, TRUE);
    
}








// ========================================================================
// initialize the population
void initialize_population(params_t *params) {

    // initialize the population from scratch

    // create initial plasmid type
    profile_t *profile = profile_new(params->beta, params->kappa, params->alpha, params);
    cell_t *cell;

    int i;

    for (i=0; i<params->psize; i++) {

	// create an empty cell
	cell = cell_new(NULL, params);
	// add the profile
	cell_add_profile(cell, profile, 1);
	// add the cell to population
	g_ptr_array_add(params->cells, cell);
	
    }

    // disown plasmid profile
    profile_release(profile);

}
