/* cell.h

   representation and functions of a cell
 
 */




#define CELL_GET_CN(cell) ((int)cell->stats->n)


// =======================================================================
//                               CELL 
// =======================================================================
typedef struct {

    // cell state variables
    double omega; 
    double omega_0;
    double domg;
    int cn; // ==> always accurate
    int age;

    double sum_alpha; // ==> always accurate

    GPtrArray *plasmids; // pointers to plasmid_t objects
    GPtrArray *mutants; // pointers to profile_t objects

    int total_extra_cn; //

    gboolean division;
    gboolean death;
    gboolean rdeath;

    nvar_t *stats; // descriptive statistics for beta, kappa, alpha

    params_t *params; // NOT OWNED BY THE CELL

} cell_t;




// create a new cell with the supplied plasmids
// the cell is reponsible for freeing the plasmid_t objects
// that are stored in cell->plasmids
cell_t *cell_new(GPtrArray *plasmids, params_t *params);


// reset the cell (to start a new cycle
void cell_reset(cell_t *cell);


// destroy the cell
void cell_free(cell_t *cell);

// add a plasmid in the cell of given profile with the supplied CN
// if the profile does not exist create a new plasmid
void cell_add_profile(cell_t *cell, profile_t *profile, int cn);

// select and return a plasmid randomly (prop. to n_i's)
plasmid_t *cell_select_plasmid(cell_t *cell);


// divide the cell and return the daughter cell
cell_t *cell_divide(cell_t *cell);


// calculate cellular growth and set the DIV/DEATH flags
// if the cell divides returns a pointer
// to the newly created daughter cell, otherwise CELL is returned
cell_t *cell_grow(cell_t *cell);


// update the cell's state (plasmid replication/conjugation)
// add conjugant plasmids into conjugation pool
// register cell in global/intra/inter stats
void cell_update(cell_t *cell, pool_t *cpool, gboolean rdeath);


// append a description of the cell to gstring
// form: (p1_descr p2_descr ...)
void cell_description(cell_t *cell, GString *st);
