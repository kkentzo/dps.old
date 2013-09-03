

// MUTATION TYPES
typedef enum {

    M_BETA,
    M_KAPPA,
    M_ALPHA,

} MUT_TYPE;



// SEGREGATION TYPES
typedef enum {

    SEG_BINOMIAL,
    SEG_PERFECT,

} SEG_TYPE;





//======================================================================
// global variable (set to true by the signal handler)
extern int global_abort;

//======================================================================
// program settings
typedef struct {

    // SIMULATION PARAMETERS
    int psize; // population size (number of cells per patch)
    int steps; // requested number of steps
    GString *log_path; // path to the log file
    GString *load_from; // file name of an HDF file to load the population from
    int seed; // the seed of the RNG
    gboolean compete; // activate competition between 2 plasmid types 
                      // requires --load_from
    int contenders; // the number of contenders (derived automatically by load_from)

    // Output parameters
    int print_every; // print info every PRINT_EVERY steps
    int fparts; // number of parts into which to split histograms
    int log_every; // how often to save stats into tables

    // MUTATIONS
    double mu; // mutation probability per cell division

    // CONJUGATION
    double pconj; // prob of conjugation (for HT_GLOBAL)

    // PLASMID PARAMETERS
    double beta; // basal replication rate
    double kappa; // sensitivity to the repressor
    double alpha; // repressor production

    double phi; // positive contribution of the plasmid trait
    double gamma; // cost coefficient for plasmid maintenance
    double gamma_alpha; // cost coefficient for inhibitor production
    double lambda; // curvature of the growth function
    double omega; // initial value for omega
    double omega_0; // basal host growth rate
    double omega_0_dev; // standard deviation of the basal host growth rate
    gboolean dilution; // whether repressors are diluted with cellular growth

    double cn_hat; // optimal CN (computed from parameters)
    int max_cn; // maximum allowed copy number (beyond that the cell dies)

    SEG_TYPE seg_type; // type of segregation during cell division

    // WHICH PARAMETERS TO MUTATE
    gboolean m_beta; // mutate beta
    gboolean m_kappa; // mutate kappa
    gboolean m_alpha; // mutate alpha

    // FORCE PARAMETER VALUES (WHEN --load_from is used)
    gboolean f_beta; // force beta
    gboolean f_kappa; // force kappa
    gboolean f_alpha; // force alpha

    // MAXIMUM PARAMETER VALUES
    double max_beta; // maximum value of beta
    double max_kappa; // maximum value of kappa
    double max_alpha; // maximum value of alpha
    int nbins; // number of bins for the 2D histograms of plasmid parameters
    
    // RANGES FOR MUTATED PARAMETERS
    double mut_rng; // the parameter (b, a) mutation range 

    // OBJECTS CREATED AT RUNTIME
    gsl_rng *rng; // random number generator object
    GArray *mwheel; // the mutation wheel

    // TIMING (set in Simulation)
    time_t start; // starting point in time of the simulation
    double duration; // set at the end of the simulation

    unsigned long long errors; // error counter

    int ndiv; // number of cell division events at the current time step
    int ndeath; // number of cell death events at the current time step
    int step; // current time step in the simulation

    // FOR USE DURING THE SIMULATION
    GPtrArray *cells; // the cellular population
    GPtrArray *dcells; // array of daughter cells
    void *logger; // the logger object
    pool_t *pool; // the plasmid profile pool


} params_t;


//======================================================================

// parse settings from command line
int parse_params(params_t *params, int argc, char **argv);

