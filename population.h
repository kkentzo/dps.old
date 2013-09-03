

// the simulation function
void run(params_t *params);



// ==================================================
//       Population Initialization and Saving
// ==================================================

// save the population of the simulation
void save_population(params_t *params, const char *fname);


// load the population from params->load_from
// and add all cells in params->cells
void load_population_from_file(params_t *params);


// initialize the population from using the default values in params
void initialize_population(params_t *params);
