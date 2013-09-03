#include <assert.h>
#include <glib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>

#include <errno.h>

#include "../utils/clib.h"
#include "../utils/nvar.h"

#include "pool.h"

#include "params.h"



#define PSIZE "psize"
#define STEPS "steps"
#define LOAD_FROM "load_from"
#define SEED "seed"
#define COMPETE "compete"

#define PRINT_EVERY "print_every"
#define FPARTS "fparts"
#define LOG_EVERY "log_every"

#define MU "mu"
#define PCONJ "pconj"
#define HT_MODEL "ht_model"

#define MUTATE "mutate"
#define FORCE "force"

#define MUT_RNG "mut_rng"

#define BETA "beta"
#define KAPPA "kappa"
#define ALPHA "alpha"

#define PHI "phi"
#define GAMMA "gamma"
#define GAMMA_ALPHA "gamma_alpha"
#define LAMBDA "lambda"
#define OMEGA "omega"
#define OMEGA_0 "omega_0"
#define OMEGA_0_DEV "omega_0_dev"
#define DILUTION "dilution"
#define SEG_TYPE "seg_type"

#define MAX_CN "max_cn"

#define MAX_BETA "max_beta"
#define MAX_KAPPA "max_kappa"
#define MAX_ALPHA "max_alpha"

#define NBINS "nbins"


//============================================================
void print_help() {

    printf("Usage: dps [options] log_path\nOptions:\n");
    printf("\nGENERAL PARAMETERS\n");
    printf("  --psize INT : set population size\n");
    printf("  --steps INT : set number of host steps\n");
    printf("  --load_from FNAME : load the population from the specified HDF file\n");
    printf("  --seed INT : set the seed of the random number generator\n");
    printf("  --print_every INT : how often (in steps) to print output \n");
    printf("  --fparts INT : number of parts into which to split histograms  \n");
    printf("  --log_every INT : frequency of writing stats to tables\n");

    printf("  -c or --compete : activate competition between 2 plasmid types\n");
    printf("                    mutations are deactivated, specify population using --load_from \n");

    printf("\nPLASMID PARAMETERS\n");
    printf("  --beta DOUBLE : basal replication rate\n");
    printf("  --kappa DOUBLE : sensitivity to the repressor\n");
    printf("  --alpha DOUBLE : repressor production\n");
    printf("  --phi DOUBLE : positive contribution of plasmid trait\n");
    printf("  --gamma DOUBLE : cost coefficient for plasmid maintenance\n");
    printf("  --gamma_alpha DOUBLE : cost coefficient for inhibitor production\n");
    printf("  --lambda DOUBLE : curvature of the host growth curve\n");
    printf("  --omega DOUBLE : initial value for omega \n");
    printf("  --omega_0 DOUBLE : mean basal host growth rate \n");
    printf("  --omega_0_dev DOUBLE : standard deviation of the basal host growth rate \n");
    printf("  --dilution [yes|no] : whether repressors are diluted with cellular growth\n");
    printf("  --seg_type [binomial|perfect] : type of plasmid segregation during cell division \n");
    printf("  --max_cn INT : maximum allowed copy number for the host\n");
    
    printf("\nMUTATIONS\n");
    printf("  --mu DOUBLE : set plasmid mutation probability (per plasmid replication)\n");
    printf("  --mut_rng DOUBLE : maximum parameter (b,k,a) mutation range \n");

    printf("  --max_beta DOUBLE : maximum allowed value of beta\n");
    printf("  --max_kappa DOUBLE : maximum allowed value of kappa\n");
    printf("  --max_alpha DOUBLE : maximum allowed value of alpha\n");
    printf("  --nbins INT : number of bins for the 2D histograms of plasmid parameters\n");

   printf("\nHORIZONTAL TRANSFER\n");
   printf("  --pconj DOUBLE : probability of conjugation (per host) for the HT_GLOBAL model\n");

    

}




//============================================================
int parse_params(params_t *params, int argc, char **argv) {
    
    int c;
    const char *optname;

    if (argc == 1) {
	printf("dps : no log path specified\nRun dps -h for usage options\n");
	return -1;
    }

    while (1) {

	int option_index = 0;
	static struct option long_options[] = {
	    {PSIZE, required_argument, 0, 0},
	    {STEPS, required_argument, 0, 0},
	    {LOAD_FROM, required_argument, 0, 0},
	    {PRINT_EVERY, required_argument, 0, 0},
	    {FPARTS, required_argument, 0, 0},
	    {LOG_EVERY, required_argument, 0, 0},
	    {SEED, required_argument, 0, 0},
	    {COMPETE, no_argument, 0, 'c'},

	    {MU, required_argument, 0, 0},
	    {PCONJ, required_argument, 0, 0},
	    {HT_MODEL, required_argument, 0, 0},
	    {MUTATE, required_argument, 0, 0},
	    {FORCE, required_argument, 0, 0},

	    {MAX_BETA, required_argument, 0, 0},
	    {MAX_KAPPA, required_argument, 0, 0},
	    {MAX_ALPHA, required_argument, 0, 0},
	    {NBINS, required_argument, 0, 0},

	    {MUT_RNG, required_argument, 0, 0},

	    {BETA, required_argument, 0, 0},
	    {KAPPA, required_argument, 0, 0},
	    {ALPHA, required_argument, 0, 0},
	    {PHI, required_argument, 0, 0},
	    {GAMMA, required_argument, 0, 0},
	    {GAMMA_ALPHA, required_argument, 0, 0},
	    {LAMBDA, required_argument, 0, 0},
	    {OMEGA, required_argument, 0, 0},
	    {OMEGA_0, required_argument, 0, 0},
	    {OMEGA_0_DEV, required_argument, 0, 0},
	    {DILUTION, required_argument, 0, 0},
	    {SEG_TYPE, required_argument, 0, 0},
	    {MAX_CN, required_argument, 0, 0},

	    {"help", no_argument, 0, 'h'},
	    // {"file", 1, 0, 0},
	    {0, 0, 0, 0}
	};

	c = getopt_long(argc, argv, "hc",
			long_options, &option_index);
	if (c == -1)
	    break;

	switch (c) {
	case 0:
	    if (optarg) {
		optname = long_options[option_index].name;

		// simulation parameters
		if (strcmp(optname, PSIZE) == 0)
		    params->psize = atoi(optarg);
	        else if (strcmp(optname, STEPS) == 0)
		    params->steps = atoi(optarg);
	        else if (strcmp(optname, LOAD_FROM) == 0)
		    params->load_from = g_string_new(optarg);
		else if (strcmp(optname, PRINT_EVERY) == 0) 
		    params->print_every = atoi(optarg);
		else if (strcmp(optname, FPARTS) == 0) 
		    params->fparts = atoi(optarg);
		else if (strcmp(optname, LOG_EVERY) == 0) 
		    params->log_every = atoi(optarg);
		
		else if (strcmp(optname, SEED) == 0) 
		    params->seed = atoi(optarg);

		// mutations
		else if (strcmp(optname, MU) == 0)
		    params->mu = atof(optarg);

		// HORIZONTAL TRANSFER
		else if (strcmp(optname, PCONJ) == 0)
		    params->pconj = atof(optarg);
		

		else if (strcmp(optname, MUTATE) == 0) {
		    if (strchr(optarg, 'b'))
			params->m_beta = TRUE;
		    if (strchr(optarg, 'k'))
			params->m_kappa = TRUE;
		    if (strchr(optarg, 'a'))
			params->m_alpha = TRUE;
		}

		else if (strcmp(optname, FORCE) == 0) {
		    if (strchr(optarg, 'b'))
			params->f_beta = TRUE;
		    if (strchr(optarg, 'k'))
			params->f_kappa = TRUE;
		    if (strchr(optarg, 'a'))
			params->f_alpha = TRUE;
		}
		
		// max param values
		else if (strcmp(optname, MAX_BETA) == 0)
		    params->max_beta = atof(optarg);
		else if (strcmp(optname, MAX_KAPPA) == 0)
		    params->max_kappa = atof(optarg);
		else if (strcmp(optname, MAX_ALPHA) == 0)
		    params->max_alpha = atof(optarg);

		else if (strcmp(optname, NBINS) == 0)
		    params->nbins = atoi(optarg);

		// maximum deviations
		else if (strcmp(optname, MUT_RNG) == 0)
		    params->mut_rng = atof(optarg);

		// plasmid parameter values
		else if (strcmp(optname, BETA) == 0)
		    params->beta = atof(optarg);
		else if (strcmp(optname, KAPPA) == 0)
		    params->kappa = atof(optarg);
		else if (strcmp(optname, ALPHA) == 0)
		    params->alpha = atof(optarg);
		else if (strcmp(optname, PHI) == 0)
		    params->phi = atof(optarg);
		else if (strcmp(optname, GAMMA) == 0)
		    params->gamma = atof(optarg);
		else if (strcmp(optname, GAMMA_ALPHA) == 0)
		    params->gamma_alpha = atof(optarg);
		else if (strcmp(optname, LAMBDA) == 0)
		    params->lambda = atof(optarg);
		else if (strcmp(optname, OMEGA) == 0)
		    params->omega = atof(optarg);
		else if (strcmp(optname, OMEGA_0) == 0)
		    params->omega_0 = atof(optarg);
		else if (strcmp(optname, OMEGA_0_DEV) == 0)
		    params->omega_0_dev = atof(optarg);
		else if (strcmp(optname, DILUTION) == 0){
		    if (strcmp(optarg, "yes") == 0)
			params->dilution = TRUE;
		    else if (strcmp(optarg, "no") == 0)
			params->dilution = FALSE;
		    else {
			printf("dilution argument not recognised : %s\nAborting\n", optarg);
			return -1;
		    }
		}
		else if (strcmp(optname, SEG_TYPE) == 0) {
		    if (strcmp(optarg, "binomial") == 0)
			params->seg_type = SEG_BINOMIAL;
		    else if (strcmp(optarg, "perfect") == 0)
			params->seg_type = SEG_PERFECT;
		    else {
			printf("seg_type argument not recognised : %s\nAborting\n", optarg);
			return -1;
		    }
		}

		else if (strcmp(optname, MAX_CN) == 0)
		    params->max_cn = atoi(optarg);
		


		printf("Setting %s=%s\n", optname, optarg);

	    }

	    break;

	case 'h':
	    print_help();
	    return 1;

	case 'c':
	    params->compete = TRUE;
	    break;
	    

	case '?':
	    return -1;

	default:
	    printf("?? getopt returned character code 0%o ??\n", c);
	}
    }

    // do we have 1 argument remaining??
    if (argc - optind != 1) {
	printf("dps: please specify a single log path (see dps -h for details)\n");
	return -1;
    }
	

    // yes we do!
    params->log_path = g_string_new(argv[optind]);

    printf("Setting log_path=%s\n", params->log_path->str);

    return 0;

}



