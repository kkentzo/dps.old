#include <string.h>
#include <assert.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <glib.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_rng.h>

#include "../utils/clib.h"
#include "../utils/nvar.h"
#include "../utils/qhdf.h"

#include "pool.h"

#include "params.h"

#include "plasmid.h"

#include "cell.h"



#include "logger.h"








// ==========================================================
// ==========================================================


// create groups -- then create tseriess
logger_t *logger_new(params_t *params) {

    logger_t *logger = malloc(sizeof(logger_t));

    // store params
    logger->params = params;

    // open file
    logger->file_id = H5Fcreate(params->log_path->str, H5F_ACC_TRUNC,
				H5P_DEFAULT, H5P_DEFAULT);

    // initialize inter- stats
    logger->inter_stats = nvar_init(P_INTER_IDX_ALL);

    // reset counters
    memset(&logger->counters, 0, sizeof(counters_t));

    // create group of dynamics
    hid_t group_id = H5Gcreate1(logger->file_id, "/dynamics", H5P_DEFAULT);
    H5Gclose(group_id);

    // ================== CREATE TABLES ===================

    // create COUNTERS table
    const char *counter_labels[] = {"n", "cn", "inf", "ptypes", "loss", "div.inf",
				    "div.all", "death", "rep", "ht", "mut"};
    hdf_table_initialize(&logger->tbl_counters, logger->file_id, "/dynamics/counters",
			 TBL_COUNTERS_LEN, Q_HDF_TYPE_INT, counter_labels);

    // create COMPETITION data??
    if (logger->params->compete) {
	// create competition group
	group_id = H5Gcreate1(logger->file_id, "/competition", H5P_DEFAULT);
	H5Gclose(group_id); 
	hdf_table_initialize(&logger->tbl_competition, logger->file_id,
			     "/competition/frequencies",
			     params->contenders, Q_HDF_TYPE_INT, NULL);
	// create contenders group
	group_id = H5Gcreate1(logger->file_id, "/competition/contenders", H5P_DEFAULT);
	H5Gclose(group_id); 

	// write all the contenders in the group
	const char *labels[] = {"beta", "kappa", "alpha"};
	double data[3];
	GHashTableIter iter;
	profile_t *profile;
	pdata_t *pdata;
	GString *name;
	hdf_table_t table;

	g_hash_table_iter_init(&iter, params->pool->profiles);

	while (g_hash_table_iter_next(&iter, (void **)&profile, (void **)&pdata)) {

	    assert(pdata->pid < params->contenders);

	    // convert profile id to a string
	    name = g_string_new_itol(pdata->pid);
	    g_string_prepend(name, "/competition/contenders/");

	    // get data into buffer
	    data[0] = profile->beta;
	    data[1] = profile->kappa;
	    data[2] = profile->alpha;

	    // create table in group contenders
	    hdf_table_initialize(&table, logger->file_id, name->str,
			     3, Q_HDF_TYPE_DOUBLE, labels);
	    hdf_table_append_record(&table, data);
	    hdf_table_finalize(&table);

	    // free string
	    g_string_free(name, TRUE);
	    
	}


    }

    nvar_t *nvar;
    
    // create INTER tables
    group_id = H5Gcreate1(logger->file_id, "/dynamics/inter", H5P_DEFAULT);
    H5Gclose(group_id);
    const char *inter_m_labels[] = {"domg", "cn", "dev", "absdev", "beta", "kappa",
				    "alpha", "bk", "rate", "nr", "ht", "death", "fitness", 
				    "tbeta", "tkappa", "talpha", "tbk"};
    hdf_table_initialize(&logger->tbl_inter_m, logger->file_id, "/dynamics/inter/M",
			 P_INTER_IDX_ALL, Q_HDF_TYPE_DOUBLE, inter_m_labels);
    hdf_table_initialize(&logger->tbl_inter_v, logger->file_id, "/dynamics/inter/V",
			 P_INTER_IDX_ALL, Q_HDF_TYPE_DOUBLE, inter_m_labels);
    // create an nvar object (just to get the covariance labels)
    nvar = nvar_init(P_INTER_IDX_ALL);
    nvar_assign_labels(nvar, inter_m_labels);
    hdf_table_initialize(&logger->tbl_inter_c, logger->file_id, "/dynamics/inter/C",
			 P_INTER_IDX_ALL*(P_INTER_IDX_ALL-1)/2, Q_HDF_TYPE_DOUBLE,
			 (const char **)nvar->plabels);
    nvar_free(nvar);

    // create INTRA tables
    group_id = H5Gcreate1(logger->file_id, "/dynamics/intra", H5P_DEFAULT);
    H5Gclose(group_id);
    const char *intra_m_labels[] = {"cn", "beta", "kappa", "alpha", "bk", "rate", 
				    "nr", "ht", "death", "fitness",
				    "tbeta", "tkappa", "talpha", "tbk"};
    hdf_table_initialize(&logger->tbl_intra_m, logger->file_id, "/dynamics/intra/M",
			 P_INTRA_IDX_ALL, Q_HDF_TYPE_DOUBLE, intra_m_labels);
    hdf_table_initialize(&logger->tbl_intra_v, logger->file_id, "/dynamics/intra/V",
			 P_INTRA_IDX_ALL, Q_HDF_TYPE_DOUBLE, intra_m_labels);
    // create an nvar object (just to get the covariance labels)
    nvar = nvar_init(P_INTRA_IDX_ALL);
    nvar_assign_labels(nvar, intra_m_labels);
    hdf_table_initialize(&logger->tbl_intra_c, logger->file_id, "/dynamics/intra/C",
			 P_INTRA_IDX_ALL*(P_INTRA_IDX_ALL-1)/2, Q_HDF_TYPE_DOUBLE,
			 (const char **)nvar->plabels);



    // create GLOBAL tables (all writing occurring within the pool)
    group_id = H5Gcreate1(logger->file_id, "/dynamics/global", H5P_DEFAULT);
    H5Gclose(group_id);
    hdf_table_initialize(&logger->tbl_global_m, logger->file_id, "/dynamics/global/M",
			 P_INTRA_IDX_ALL, Q_HDF_TYPE_DOUBLE, intra_m_labels);
    hdf_table_initialize(&logger->tbl_global_v, logger->file_id, "/dynamics/global/V",
			 P_INTRA_IDX_ALL, Q_HDF_TYPE_DOUBLE, intra_m_labels);
    hdf_table_initialize(&logger->tbl_global_c, logger->file_id, "/dynamics/global/C",
			 P_INTRA_IDX_ALL*(P_INTRA_IDX_ALL-1)/2, Q_HDF_TYPE_DOUBLE,
			 (const char **)nvar->plabels);
    // free the nvar object
    nvar_free(nvar);
    
    // reset the intra-buffers
    BFILL(logger->MM, P_INTRA_IDX_ALL, 0);
    BFILL(logger->VV, P_INTRA_IDX_ALL, 0);
    BFILL(logger->CC, P_INTRA_IDX_ALL * (P_INTRA_IDX_ALL-1) / 2, 0);

    // === CREATE HISTOGRAMS ===
    // create group
    group_id = H5Gcreate1(logger->file_id, "/histograms", H5P_DEFAULT);
    H5Gclose(group_id);

    // initialize CN histogram
    double xlim[] = {-0.5, params->max_cn - 0.5};
    hdf_histogram_initialize(&logger->hist_cn, logger->file_id, "/histograms/cn",
			     params->fparts, params->max_cn, 0,
			     xlim, NULL, FALSE);

    // initialize AGE histogram
    xlim[0] = 0;
    xlim[1] = 50;
    hdf_histogram_initialize(&logger->hist_age, logger->file_id, "/histograms/age",
			     params->fparts, 50, 0,
			     xlim, NULL, FALSE);

    // initialize BK histogram
    xlim[1] = logger->params->max_beta;
    double ylim[] = {0, logger->params->max_kappa};
    hdf_histogram_initialize(&logger->hist_beta_kappa, logger->file_id, "/histograms/bk",
			     params->fparts, logger->params->nbins, logger->params->nbins,
			     xlim, ylim, TRUE);

    // initialize BA histogram
    xlim[0] = 0; xlim[1] = logger->params->max_beta;
    ylim[0] = 0; ylim[1] = logger->params->max_alpha;
    hdf_histogram_initialize(&logger->hist_beta_alpha, logger->file_id, "/histograms/ba",
			     params->fparts, logger->params->nbins, logger->params->nbins,
			     xlim, ylim, TRUE);


    // initialize KA histogram
    xlim[0] = 0; xlim[1] = logger->params->max_kappa;
    ylim[0] = 0; ylim[1] = logger->params->max_alpha;
    hdf_histogram_initialize(&logger->hist_kappa_alpha, logger->file_id, "/histograms/ka",
			     params->fparts, logger->params->nbins, logger->params->nbins,
			     xlim, ylim, TRUE);

    
    return logger;
    
}




void logger_log(logger_t *logger, int step) {

    // write ptypes
    logger->counters.ptypes = logger->params->pool->size;

    //assert(logger->counters.cn == logger->inter_stats->n);
    
    // print info on screen??
    if (logger->params->print_every && (step % logger->params->print_every == 0))
    	printf("step %5d | %4d/%d (%d) || b=%.3f | k=%.3f | a=%.3f || cn=%.3f\n",
    	       step, logger->counters.inf, logger->counters.n, 
	       logger->counters.ptypes,
	       logger->beta,
	       logger->kappa,
	       logger->alpha,
	       1. * logger->counters.cn / logger->counters.inf);


    if (step % logger->params->log_every == 0) {

	// ================ INTER ==================
	// write the inter-stats to buffers
	nvar_get_statistics(logger->inter_stats, logger->M, logger->V, logger->C,
			    FALSE, 1);
	// write the buffers to the tables
	hdf_table_append_record(&logger->tbl_inter_m, logger->M);
	hdf_table_append_record(&logger->tbl_inter_v, logger->V);
	hdf_table_append_record(&logger->tbl_inter_c, logger->C);
	// reset the stats
	nvar_reset(logger->inter_stats);
	// ================ INTRA ==================
	// calculate means
	int i;
	for (i=0; i<P_INTRA_IDX_ALL; i++) {
	    
	    logger->MM[i] /= logger->counters.cn;
	    logger->VV[i] /= logger->counters.cn;

	}
	for (i=0; i<P_INTRA_IDX_ALL*(P_INTRA_IDX_ALL-1)/2; i++)
	    logger->CC[i] /= logger->counters.cn;

	// write the buffers to the tables
	hdf_table_append_record(&logger->tbl_intra_m, logger->MM);
	hdf_table_append_record(&logger->tbl_intra_v, logger->VV);
	hdf_table_append_record(&logger->tbl_intra_c, logger->CC);
	// reset the buffers
	BFILL(logger->MM, P_INTRA_IDX_ALL, 0);
	BFILL(logger->VV, P_INTRA_IDX_ALL, 0);
	BFILL(logger->CC, P_INTRA_IDX_ALL * (P_INTRA_IDX_ALL-1) / 2, 0);

	// write counters to table
	hdf_table_append_record(&logger->tbl_counters, &logger->counters);
	// reset counters
	memset(&logger->counters, 0, sizeof(counters_t));
	
    }

    // write and reset histograms??
    int save_every = logger->params->steps / logger->params->fparts;
    if (step > 0 && save_every > 0 && step % save_every == 0) {
	hdf_histogram_write(&logger->hist_cn);
	hdf_histogram_write(&logger->hist_age);
	hdf_histogram_write(&logger->hist_beta_kappa);
	hdf_histogram_write(&logger->hist_beta_alpha);
	hdf_histogram_write(&logger->hist_kappa_alpha);
    }

}






// record fitness, infection and register copy number in histogram
void logger_register_cell(logger_t *logger, cell_t *cell) {

    // count cell towards population size
    logger->counters.n += 1;
    logger->counters.cn += cell->cn;

    // register copy number
    hdf_histogram_add(&logger->hist_cn, cell->cn);

    if (cell->cn > 0) {

	++ logger->counters.inf;

	double host[] = {cell->domg,
			 cell->cn,
			 cell->cn - logger->params->cn_hat,
			 abs(cell->cn - logger->params->cn_hat),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_BETA),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_KAPPA),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_ALPHA),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_BK),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_RATE),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_N_R),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_HT),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_DEATH),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_FITNESS),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_TBETA),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_TKAPPA),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_TALPHA),
			 nvar_calc_mean(cell->stats, P_INTRA_IDX_TBK)};

	// register cell in inter- stats
	nvar_add(logger->inter_stats, host, cell->cn);

	// update the intra-stats buffers
	// add the means, variances and covariances weighted by cell's CN
	nvar_get_statistics(cell->stats, logger->MM, logger->VV, logger->CC,
			    TRUE, cell->cn);
    }
}







// register a cell division event
void logger_register_division(logger_t *logger, cell_t *cell) {

    // is the cell infected??
    if (cell->cn > 0) {

	// register division of a plasmid-infected cell
	logger_register_events(logger, E_DIV_INF, 1);
	
    }

    // register division
    logger_register_events(logger, E_DIV_ALL, 1);

    // register cell's age
    hdf_histogram_add(&logger->hist_age, cell->age);
    
}








// register a number of events for event type ETYPE
void logger_register_events(logger_t *logger, EVENT_T etype, int events) {

    switch (etype)
	{
	case E_LOSS :
	    logger->counters.loss += events; return;
	case E_DIV_INF :
	    logger->counters.div_inf += events; return;
	case E_DIV_ALL :
	    logger->counters.div_all += events; return;
	case E_DEATH :
	    logger->counters.death += events; return;
	case E_REP :
	    logger->counters.rep += events; return;
	case E_HT :
	    logger->counters.ht += events; return;
	case E_MUT :
	    logger->counters.mut += events; return;
	    
	    
	}

}





void logger_write_params(logger_t *logger) {

    // create settings group
    hid_t group_id = H5Gcreate1(logger->file_id, "/settings", H5P_DEFAULT);
    H5Gclose(group_id);    

    // create dataspace
    hsize_t dim = 1;

    H5LTmake_dataset(logger->file_id, "/settings/psize", 1, &dim,
		     H5T_NATIVE_INT, &logger->params->psize);

    H5LTmake_dataset(logger->file_id, "/settings/steps", 1, &dim,
		     H5T_NATIVE_INT, &logger->params->steps);

    H5LTmake_dataset(logger->file_id, "/settings/fparts", 1, &dim,
		     H5T_NATIVE_INT, &logger->params->fparts);
    
    H5LTmake_dataset(logger->file_id, "/settings/seed", 1, &dim,
    		     H5T_NATIVE_INT, &logger->params->seed);
    
    H5LTmake_dataset(logger->file_id, "/settings/seg_type", 1, &dim,
    		     H5T_NATIVE_INT, &logger->params->seg_type);

    H5LTmake_dataset(logger->file_id, "/settings/mu", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->mu);

    H5LTmake_dataset(logger->file_id, "/settings/pconj", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->pconj);

    H5LTmake_dataset(logger->file_id, "/settings/beta", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->beta);

    H5LTmake_dataset(logger->file_id, "/settings/kappa", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->kappa);
    
    H5LTmake_dataset(logger->file_id, "/settings/alpha", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->alpha);

    H5LTmake_dataset(logger->file_id, "/settings/max_beta", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->max_beta);

    H5LTmake_dataset(logger->file_id, "/settings/max_kappa", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->max_kappa);
    
    H5LTmake_dataset(logger->file_id, "/settings/max_alpha", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->max_alpha);

    H5LTmake_dataset(logger->file_id, "/settings/max_cn", 1, &dim,
    		     H5T_NATIVE_INT, &logger->params->max_cn);
    
    H5LTmake_dataset(logger->file_id, "/settings/phi", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->phi);

    H5LTmake_dataset(logger->file_id, "/settings/gamma", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->gamma);

    H5LTmake_dataset(logger->file_id, "/settings/gamma.alpha", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->gamma_alpha);
    
    H5LTmake_dataset(logger->file_id, "/settings/lambda", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->lambda);

    H5LTmake_dataset(logger->file_id, "/settings/omega.0", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->omega_0);

    H5LTmake_dataset(logger->file_id, "/settings/dilution", 1, &dim,
    		     H5T_NATIVE_INT, &logger->params->dilution);

    H5LTmake_dataset(logger->file_id, "/settings/mut.rng", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->mut_rng);

    H5LTmake_dataset(logger->file_id, "/settings/duration", 1, &dim,
    		     H5T_NATIVE_DOUBLE, &logger->params->duration);
    
    

}



// save the population
void logger_save_population(logger_t *logger, GString *profiles, GString *hosts) {

    // create population group
    hid_t group_id;
    group_id = H5Gcreate1(logger->file_id, "/population", H5P_DEFAULT);

    H5LTmake_dataset_string(group_id, "/population/profiles", profiles->str);
    H5LTmake_dataset_string(group_id, "/population/hosts", hosts->str);

    H5Gclose(group_id);
    
}




void logger_close(logger_t *logger) {

    // close the histograms

    hdf_histogram_finalize(&logger->hist_cn);
    hdf_histogram_finalize(&logger->hist_age);
    hdf_histogram_finalize(&logger->hist_beta_kappa);
    hdf_histogram_finalize(&logger->hist_beta_alpha);
    hdf_histogram_finalize(&logger->hist_kappa_alpha);

    // close the tables
    hdf_table_finalize(&logger->tbl_counters);

    if (logger->params->compete) 
	hdf_table_finalize(&logger->tbl_competition);


    hdf_table_finalize(&logger->tbl_global_m);
    hdf_table_finalize(&logger->tbl_global_v);
    hdf_table_finalize(&logger->tbl_global_c);    

    hdf_table_finalize(&logger->tbl_inter_m);
    hdf_table_finalize(&logger->tbl_inter_v);
    hdf_table_finalize(&logger->tbl_inter_c);    

    hdf_table_finalize(&logger->tbl_intra_m);
    hdf_table_finalize(&logger->tbl_intra_v);
    hdf_table_finalize(&logger->tbl_intra_c);    
    
    // Write the parameter values
    logger_write_params(logger);

    // close stats
    nvar_free(logger->inter_stats);

    // close file
    H5Fclose(logger->file_id);

    // free the struct
    free(logger);
    
}
