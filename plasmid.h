/*

  plasmid.h

 */



// ==================================================================================
//                          A PARASITE PROFILE
//    characterized by a particular configuration of beta and alpha
// ==================================================================================


typedef struct {

    // profile
    double beta; // basal replication rate
    double kappa; // sensitivity to the repressor
    double alpha; // repressor production

    // object's reference count
    // represents how many hosts the profile is spread across
    int _ref_count;

    // ref to parent profile
    struct profile_t *parent;

    // mark the profile creation sim. step
    int step;

    // ref to params
    params_t *params;

} profile_t;



// create a new plasmid replication profile
profile_t * profile_new(double beta, double kappa, double alpha, params_t *params);


// mutate the profile and return a new one
profile_t * profile_new_m(profile_t *profile);

// set the profile's parent (used in profile_new_m())
void profile_set_parent(profile_t *profile, profile_t *parent);


// increase the object's reference count
profile_t *profile_retain(profile_t *profile);


// reduce the object's reference count
// if it goes to zero the memory is freed as well
profile_t *profile_release(profile_t *profile);


// return the global transmission rate of the profile
//double profile_calc_transmission_rate(profile_t *profile);


// increase/decrease the profile's copy number in the population
void profile_inc(profile_t *profile, int step);
void profile_dec(profile_t *profile, int step);


// return the total cn of the profile
int profile_get_cn(profile_t *profile);



// append the profile's description in the supplied string
// form: (mem_add beta alpha ref_count)
void profile_description(profile_t *profile, GString *st);




// ========================================================================
//                      A PLASMID IN THE CELL
//             characterized by a profile and a copy number
// ========================================================================


typedef struct {

    profile_t *profile;

    int cn;

} plasmid_t;



// create a new plasmid with the supplied profile and copy number
// ** retains the profile **
plasmid_t * plasmid_new(profile_t *profile, int cn);


// increase/decrease the plasmid's copy number
void plasmid_inc(plasmid_t *plasmid, int step);
void plasmid_dec(plasmid_t *plasmid, int step);


// free the plasmid 
// ** release the profile **
void plasmid_free(plasmid_t *plasmid);


// append the plasmid's description to the gstring
// form: (cn profile_mem_address)
void plasmid_description(plasmid_t *plasmid, GString *st);

