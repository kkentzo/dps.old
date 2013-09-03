




// the indices of means (and vars) of dynamics/[global|inter|intra]


// used for INTRA and GLOBAL
typedef enum {

    P_INTRA_IDX_CN,
    P_INTRA_IDX_BETA,
    P_INTRA_IDX_KAPPA,
    P_INTRA_IDX_ALPHA,
    P_INTRA_IDX_BK,
    P_INTRA_IDX_RATE,    
    P_INTRA_IDX_N_R,
    P_INTRA_IDX_HT,
    P_INTRA_IDX_DEATH,
    P_INTRA_IDX_FITNESS,
    P_INTRA_IDX_TBETA,
    P_INTRA_IDX_TKAPPA,
    P_INTRA_IDX_TALPHA,
    P_INTRA_IDX_TBK,
    
    P_INTRA_IDX_ALL,

} P_INTRA_IDX;


// used for INTER
typedef enum {

    P_INTER_IDX_DOMG,
    P_INTER_IDX_CN,
    P_INTER_IDX_CN_DEV,
    P_INTER_IDX_CN_ABS_DEV,
    P_INTER_IDX_BETA,
    P_INTER_IDX_KAPPA,
    P_INTER_IDX_ALPHA,
    P_INTER_IDX_BK,
    P_INTER_IDX_RATE,
    P_INTER_IDX_N_R,
    P_INTER_IDX_HT,
    P_INTER_IDX_DEATH,
    P_INTER_IDX_FITNESS,
    P_INTER_IDX_TBETA,
    P_INTER_IDX_TKAPPA,
    P_INTER_IDX_TALPHA,
    P_INTER_IDX_TBK,
    
    P_INTER_IDX_ALL,

} P_INTER_IDX;







// EVENT TYPES
typedef enum {
    
    E_LOSS,
    E_DIV_INF,
    E_DIV_ALL,
    E_DEATH,
    E_REP,
    E_HT,
    E_MUT,

} EVENT_T;









// COUNTER STRUCT
typedef struct {

    int n;
    int cn;
    int inf;
    int ptypes;
    int loss;
    int div_inf;
    int div_all;
    int death;
    int rep;
    int ht;
    int mut;

} counters_t;


#define TBL_COUNTERS_LEN 11


// ==========================================================
// ==========================================================

typedef struct {

    params_t *params;

    hid_t file_id;

    counters_t counters;

    // plasmid parameter means : set directly in pool_update()
    double beta;
    double kappa;
    double alpha;

    // tables
    hdf_table_t tbl_counters;

    hdf_table_t tbl_competition;

    hdf_table_t tbl_global_m;
    hdf_table_t tbl_global_v;
    hdf_table_t tbl_global_c;

    hdf_table_t tbl_inter_m;
    hdf_table_t tbl_inter_v;
    hdf_table_t tbl_inter_c;

    hdf_table_t tbl_intra_m;
    hdf_table_t tbl_intra_v;
    hdf_table_t tbl_intra_c;

    // buffers for transferring inter-cellular statistics to inter- tables
    nvar_t *inter_stats;
    double M[P_INTER_IDX_ALL];
    double V[P_INTER_IDX_ALL];
    double C[P_INTER_IDX_ALL*(P_INTER_IDX_ALL-1)/2];

    // buffers for adding intra- stats
    double MM[P_INTRA_IDX_ALL];
    double VV[P_INTRA_IDX_ALL];
    double CC[P_INTRA_IDX_ALL*(P_INTRA_IDX_ALL-1)/2];

    // datasets for histograms
    // for each histogram remember to :
    // 1. create the group in logger_new()
    // 2. call hdf_histogram_initialize() in logger_new()
    // 3. call hdf_histogram_inform() in logger_log()
    // 4. call hdf_histogram_finalize() in logger_close()
    hdf_histogram_t hist_cn;
    hdf_histogram_t hist_age;
    hdf_histogram_t hist_beta_kappa;
    hdf_histogram_t hist_beta_alpha;
    hdf_histogram_t hist_kappa_alpha;


} logger_t;



// create groups -- then create tseries
logger_t *logger_new(params_t *params);


// update everything and move to the next step
void logger_log(logger_t *logger, int step);


// record infection and register copy number in histogram
void logger_register_cell(logger_t *logger, cell_t *cell);

// register a cell division event
void logger_register_division(logger_t *logger, cell_t *cell);

// register a number of events for event type ETYPE
void logger_register_events(logger_t *logger, EVENT_T etype, int events);

// save the population
void logger_save_population(logger_t *logger, GString *profiles, GString *hosts);

void logger_close(logger_t *logger);
