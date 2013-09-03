
// ==========================================
// THE GLOBAL POOL OF UNIQUE PLASMID PROFILES
// ==========================================



typedef struct {

    long pid; // the profile's pid

    int cn; // the profile's total copy number

    // the sums of the mutant profile values
    double mbeta_sum;
    double mkappa_sum;
    double malpha_sum;
    double mbk_sum;
    // FITNESS_SUM is the sum of the profiles different fitnesses
    // across cells (weighted by copy number)
    // It represents the total representation of the type
    // at the next time step (i.e. it is the number of offsrping)
    double fitness_sum; 

    int pcn; // the copy number of the profile at the current time step
    int mutants; // the number of mutation events

} pdata_t;




typedef struct {

    GHashTable *profiles;

    long pid; // the id of the next plasmid profile to be created

    nvar_t *stats; // pool statistics

    // the total copy number
    int cn;

    // the number of unique plasmid profiles
    int size;

    // buffers for means, vars and covs
    double *M;
    double *V;
    double *C;


} pool_t;




// initialize the pool
pool_t *pool_new();

// free the pool
void pool_free(pool_t *pool);


// empty the pool
void pool_empty(pool_t *pool);


// Add a profile to the pool (with CN=0)
void pool_add_profile(pool_t *pool, void *profile);

// remove a profile 
void pool_remove_profile(pool_t *pool, void *profile);

// increase/decrease the copy number of a profile
// ** the profile is added to the pool if it's not there already **
void pool_change_profile_count(pool_t *pool, void *profile, int step);


// return the CN of a profile
int pool_get_profile_cn(pool_t *pool, void *profile);

// register intra-cellular stats
void pool_register_cell_stats(pool_t *pool, nvar_t *stats);

// register the profile from within cell_update() (for computing transmission biases)
void pool_register_profile(pool_t *pool, void *profile, int cn,
			   double fitness, int mutants);

// register a mutation event (for computing transmission biases)
void pool_register_mutation(pool_t *pool, void *wt_profile, void *mut_profile);

void pool_description(pool_t *pool, GString *st);

// reset losses for all profiles and update joint distributions in logger
void pool_update(pool_t *pool, int step, void *_params);

// distribute the conjugants among recipients
void pool_distribute(pool_t *pool, GPtrArray *recipients, void *_params);
