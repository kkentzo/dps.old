#include <hdf5.h>
#include <glib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <sys/stat.h> // for mkdir()

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>


#include "../utils/clib.h"
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



//============================================================
// default settings
params_t params = {

    // SIMULATION PARAMETERS
    1000, // psize
    50000, // steps
    NULL, // log_path
    NULL, // load_from
    0, // seed
    FALSE, // compete
    1, // contenders

    // Output parameters
    5000, // print_every
    20, // fparts
    1, // log_every

    5e-3, // mu (mutation)

    0, // pconj

    // PLASMID PARAMETERS
    0.05, // beta
    0, // kappa
    0, // alpha

    0.1, // phi
    3e-3, // gamma
    0, // gamma_alpha
    3, // lambda
    1, // omega
    0.05, // omega_0
    0, // omega_0_dev
    TRUE, // dilution

    0, // cn_hat
    70, // max_cn    

    SEG_BINOMIAL, // seg_type

    // WHICH PARAMETERS TO MUTATE
    FALSE, // m_beta
    FALSE, // m_kappa
    FALSE, // m_alpha

    // FORCE PARAMETER VALUES (WHEN --load_from is used)
    FALSE, // f_beta
    FALSE, // f_kappa
    FALSE, // f_alpha
    
    // MAXIMUM PARAMETER VALUES
    1, // max_beta
    1, // max_kappa
    1, // max_alpha
    20, // nbins

    // MAXIMUM DEVIATIONS FOR MUTATED PARAMETERS
    0.05, // mut_rng

    // OBJECTS CREATED AT RUNTIME
    NULL, // rng
    NULL, // mwheel

    // TIMING (set in Simulation)
    0, // start
    0, // duration

    // counters
    0, // error counter

    0, // ndiv
    0, // ndeath
    0, // step

    NULL, // cells
    NULL, // dcells
    NULL, // logger
    NULL, // profiles

};






//======================================================
// create all the objects that are needed by the process
void prepare(params_t *params) {


    // === initialize RNG ===
    // init default generator
    params->rng = gsl_rng_alloc(gsl_rng_default);
    // seed the generator
    if (params->seed == 0)
	// draw a seed from /dev/urandom
	params->seed = rand();
    printf("Seed = %d\n", params->seed);
    // set the seed
    gsl_rng_set(params->rng, params->seed);

    // initialize the mutation wheel
    int var;
    params->mwheel = g_array_new(FALSE, FALSE, sizeof(int));
    if (params->m_beta) {
	var = M_BETA;
	g_array_append_val(params->mwheel, var);
    }
    if (params->m_kappa) {
	var = M_KAPPA;
	g_array_append_val(params->mwheel, var);
    }
    if (params->m_alpha) {
	var = M_ALPHA;
	g_array_append_val(params->mwheel, var);
    }

    // compute the optimal copy number
    params->cn_hat = sqrt(params->phi * params->lambda / params->gamma) 
	- params->lambda;

    printf("Optimal CN is %.2f\n", params->cn_hat);


    // initialize the population
    params->cells = g_ptr_array_new();
    // initialize array of daughter cells
    params->dcells = g_ptr_array_new();
    // initialize the profile pool
    params->pool = pool_new();


    // load the population from a file??
    if (params->load_from)
	load_population_from_file(params);
    else
	initialize_population(params);

    // do we have a competition??
    if (params->compete) {

	// set the number of contenders based on the pool's contents
	params->contenders = params->pool->size;

	printf("COMPETITION MODE : %d contenders found\n", params->contenders);

	// OK, we do -- now cancel all mutations
	g_array_set_size(params->mwheel, 0);
	
    }

    // adjust steps so that it is a multiple of fparts
    int steps = params->steps;
    while (params->steps % params->fparts > 0)
	++ params->steps;
    if (steps != params->steps)
	printf("Adjusted simulation steps so that it is a multiple of FPARTS (old=%d, new=%d)\n",
	       steps, params->steps);

    // initialize logger
    params->logger = logger_new(params);


}




//=====================================
// save some stuff and free all objects
void finalize(params_t *params) {

    // free objects
    gsl_rng_free(params->rng);
    g_array_free(params->mwheel, TRUE);

    // close the stats object
    logger_close(params->logger);

    // free all cells in the population
    int i;
    for (i=0; i<params->cells->len; i++) 
	cell_free(params->cells->pdata[i]);
    
    // free cell array
    g_ptr_array_free(params->cells, TRUE);
    // free the array of daughter cells
    g_ptr_array_free(params->dcells, TRUE);
    // free the profile pool
    pool_free(params->pool);

    // ** the log_path will be freed in the end of main **

}




// ============================================================
// THIS WILL BE SET TO TRUE BY THE SIGNAL HANDLER

int global_abort = 0;

// C-c handler
void handler(int signal) {

    global_abort = 1;

    printf("\nReceived signal %d\n", signal);

}





//============================================

int main(int argc, char **argv) {


    // setup gsl environment
    gsl_rng_env_setup();
 
    // register C-c handler
    signal(SIGINT, handler);
   

    // parse cmd line args (overrides default params)
    int status = parse_params(&params, argc, argv);
    if (status) 
    	return -1;

    // initialize the necessary objects in params
    prepare(&params);

    // get and print time
    time_t rawtime;
    time(&rawtime);
    printf("Starting simulation :: %s", ctime(&rawtime));
    
    // run the simulation
    run(&params);

    // print duration
    GString *dest = g_string_new("");
    sec_to_str(params.duration, dest);
    printf("Process finished (duration : %s)\n", dest->str);
    g_string_free(dest, TRUE);

    // release all objects in params
    finalize(&params);

    // free remaining strings
    if (params.load_from)
	g_string_free(params.load_from, TRUE);
    
    g_string_free(params.log_path, TRUE);

    return 0;

}
