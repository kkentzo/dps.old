#include <stdlib.h>
#include <string.h>
#include <glib.h>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include <hdf5.h>
#include <hdf5_hl.h>
#include "qhdf.h"
#include "clib.h"






//======================================================================================
//======================================================================================


// initialize a table (struct already exists)
// FIELD_NAMES can be NULL, in which case names are auto-generated as letters of
// the alphabet; the case were nfields >= LETTERS_SIZE is also handled (numbers are
// appended to letters accordingly)
void hdf_table_initialize(hdf_table_t *table, hid_t file_id, const char *dataset_name, 
			  size_t nfields, Q_HDF_TYPE table_type, const char **field_names)
{

    // store pointers to constant strings
    table->file_id = file_id;
    table->dataset_name = dataset_name;

    table->_chunk_size_ = TABLE_DEFAULT_CHUNK_SIZE;

    // determine field size and hdf type
    size_t field_size;
    hid_t field_hdf_type;

    switch (table_type)
	{
	case Q_HDF_TYPE_INT :
	    field_size = sizeof(int);
	    field_hdf_type = H5T_NATIVE_INT;
	    break;
	case Q_HDF_TYPE_DOUBLE :
	    field_size = sizeof(double);
	    field_hdf_type = H5T_NATIVE_DOUBLE;
	    break;
	}

    table->nfields = nfields;
    table->record_size = nfields * field_size;

    int i;
    // create the types, offsets, the sizes and store the labels
    table->field_types = malloc(nfields * sizeof(hid_t));
    table->field_sizes = malloc(nfields * sizeof(size_t));
    table->field_offsets = malloc(nfields * sizeof(size_t));
    table->field_names = malloc(nfields * sizeof(char *));

    // should we automaticall generate field names?
    table->generated_names = NULL;
    if (field_names == NULL) {

	GString *name;

	// create list of strings
	table->generated_names = g_ptr_array_new();

	for (i=0; i<nfields; i++) {
	    // create string
	    name = g_string_new_itol(i);
	    // add it to array
	    g_ptr_array_add(table->generated_names, name);
	    // add the CString pointer to field_names
	    table->field_names[i] = name->str;
	}
	
    }
    
    
    for (i=0; i<nfields; i++) {
	table->field_types[i] = field_hdf_type;
	table->field_sizes[i] = field_size;
	table->field_offsets[i] = i * field_size;
	if (field_names)
	    table->field_names[i] = field_names[i];
    }

    // create the table in the file (no data are actually written)
    H5TBmake_table("Table", table->file_id, table->dataset_name,
		   table->nfields, 0, table->record_size,
		   table->field_names, table->field_offsets,
		   table->field_types, table->_chunk_size_,
		   NULL, TRUE, NULL);

    // create the record buffer
    table->buf_size = TABLE_DEFAULT_CHUNK_SIZE;
    table->buf = calloc(table->buf_size, table->record_size);
    table->buf_idx = 0;
    
}







void hdf_table_append_record(hdf_table_t *table, void *data) {

    if (table->buf_idx == table->buf_size)
	hdf_table_flush(table);

    // copy the data into the buffer
    memcpy(table->buf + table->record_size * table->buf_idx, data, table->record_size);
    table->buf_idx ++;

}






void hdf_table_flush(hdf_table_t *table) {

    H5TBappend_records(table->file_id, table->dataset_name,
		       table->buf_idx, table->record_size,
		       table->field_offsets, table->field_sizes,
		       table->buf);

    // reset the buffer index
    table->buf_idx = 0;
}







// attach (write) an attribute to the table
/* void hdf_table_attach_attribute(hdf_table_t *table) { */

/*     hid_t dataset_id, dataspace_id, attribute_id; */
/*     hsize_t dims = 3; */
/*     double data[dims]; // beta, kappa, alpha values for each type */
/*     data[0] = 1; */
/*     data[1] = 2; */
/*     data[2] = 3; */

/*     /\* switch (table_type) *\/ */
/*     /\* 	{ *\/ */
/*     /\* 	case Q_HDF_TYPE_INT : *\/ */
/*     /\* 	    field_size = sizeof(int); *\/ */
/*     /\* 	    field_hdf_type = H5T_NATIVE_INT; *\/ */
/*     /\* 	    break; *\/ */
/*     /\* 	case Q_HDF_TYPE_DOUBLE : *\/ */
/*     /\* 	    field_size = sizeof(double); *\/ */
/*     /\* 	    field_hdf_type = H5T_NATIVE_DOUBLE; *\/ */
/*     /\* 	    break; *\/ */
/*     /\* 	} *\/ */
    

/*     // open dataset */
/*     dataset_id = H5Dopen2(table->file_id, table->dataset_name, H5P_DEFAULT); */
/*     // define dataspace */
/*     dataspace_id = H5Screate_simple(1, &dims, NULL); */
/*     // create attribute */
/*     attribute_id = H5Acreate2(dataset_id, "attr", H5T_NATIVE_DOUBLE, dataspace_id, */
/* 			      H5P_DEFAULT, H5P_DEFAULT); */
/*     // write attribute data */
/*     H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, data); */

/*     // close objects */
/*     H5Aclose(attribute_id); */
/*     H5Sclose(dataspace_id); */
/*     H5Dclose(dataset_id); */
    
/* } */






void hdf_table_finalize(hdf_table_t *table) {

    // flush the buffer
    if (table->buf_idx > 0)
	hdf_table_flush(table);

    // free the resources
    free(table->field_types);
    free(table->field_sizes);
    free(table->field_offsets);
    free(table->field_names);

    // free generated names??
    if (table->generated_names != NULL) {

	// free all the strings
	int i;
	for (i=0; i<table->generated_names->len; i++)
	    g_string_free(g_ptr_array_index(table->generated_names, i), TRUE);

	// free the pointer array
	g_ptr_array_free(table->generated_names, TRUE);
	
    }

    free(table->buf);    
    
}








// ===========================================================================
//        READ A TABLE FROM A FILE AND RETURN IT IN MATRIX FORM
// ===========================================================================
// loads the requested table from open file FILE_ID
// OUT: writes the number of records in NRECORDS
// returns a newly allocated buffer with the table's contents
// remember to FREE the buffer
void *hdf_table_read(hid_t file_id, const char *table_name, hsize_t *nrecords) {

    // read profiles table
    hsize_t nfields;
    H5TBget_table_info(file_id, table_name, &nfields, nrecords);

    printf("Loading table %s : %d records (fields=%d)\n",
	   table_name, (int)*nrecords, (int)nfields);

    // get field info
    //char *field_names[nfields];
    /* int i; */
    /* for (i=0; i<nfields; i++) */
    /* 	field_names[i] = malloc(max_name_size * sizeof(char)); */
    size_t field_sizes[nfields];
    size_t field_offsets[nfields];
    size_t type_size;
    H5TBget_field_info(file_id, table_name, NULL, field_sizes,
		       field_offsets, &type_size);

    // create table in memory
    void *buffer = malloc((*nrecords) * type_size);

    H5TBread_table(file_id, table_name, type_size, field_offsets,
		   field_sizes, buffer);

    return buffer;
    
}






//======================================================================================
//======================================================================================



void hdf_dataset_initialize(hdf_dataset_t *dataset, hid_t file_id, const char *name,
			    hsize_t ndim, hsize_t *dims, Q_HDF_TYPE dataset_type)
{

    // store file ID
    dataset->file_id = file_id;

    // store dataset type
    switch (dataset_type)
	{
	case Q_HDF_TYPE_INT :
	    dataset->dataset_type = H5T_NATIVE_INT;
	    break;
	case Q_HDF_TYPE_DOUBLE :
	    dataset->dataset_type = H5T_NATIVE_DOUBLE;
	    break;
	}

    // store number of dimensions
    dataset->ndim = ndim;

    // store dimensionality
    dataset->dims = malloc(dataset->ndim * sizeof(hsize_t));
    memmove(dataset->dims, dims, dataset->ndim * sizeof(hsize_t));


    // create dataspace
    hid_t dataspace_id = H5Screate_simple(dataset->ndim, dataset->dims, NULL);

    hid_t plist_id = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, dataset->ndim, dataset->dims);
    H5Pset_deflate(plist_id, 6);
    // create the dataset
    dataset->dataset_id = H5Dcreate2(dataset->file_id, name, dataset->dataset_type,
				     dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);

    H5Pclose(plist_id);

    // close the dataspace
    H5Sclose(dataspace_id);
    
}





// close the dataset
void hdf_dataset_finalize(hdf_dataset_t *dataset) {

    // close the dataset
    H5Dclose(dataset->dataset_id);

    // free the dimensionality buffer
    free(dataset->dims);

}





// TR_SIZE is a buffer of size dataset->ndim which contains
// the truncation sizes along each dimension
void hdf_dataset_truncate(hdf_dataset_t *dataset, hsize_t *tr_size) {

    if (tr_size) 
	H5Dset_extent(dataset->dataset_id, tr_size);

}




// write the DATA into the dataset -- 
// pass NULL to all to write to the whole dataset
// ** DATA are written along the rightmost dimension first (i.e. byrow=T) **
void hdf_dataset_write(hdf_dataset_t *dataset, void *data, hsize_t *offsets, hsize_t *counts)
{

    hid_t dataspace_id, memspace_id;

    hsize_t *dims = ((offsets && counts) ? counts : dataset->dims);

    // get memspace and dataspace
    memspace_id = H5Screate_simple(dataset->ndim, dims, NULL);
    dataspace_id = H5Dget_space(dataset->dataset_id);

    // select hyperslab within the dataset??
    if (offsets && counts)
	H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offsets,
			    NULL, counts, NULL);

    // write the data
    H5Dwrite(dataset->dataset_id, dataset->dataset_type, memspace_id,
	     dataspace_id, H5P_DEFAULT, data);

    // close spaces
    H5Sclose(dataspace_id);
    H5Sclose(memspace_id);
    
}




// simplified version of dataset_write()
void hdf_dataset_write_all(hdf_dataset_t *dataset, void *data) {

    hdf_dataset_write(dataset, data, NULL, NULL);
    
}





//======================================================================================
//======================================================================================




// create the objects

void hdf_histogram_initialize(hdf_histogram_t *hist, hid_t file_id, const char *path, 
			      int n_snapshots, hsize_t nbins1, hsize_t nbins2,
			      double *xlim, double *ylim,
			      int adjust_right_bin)
{

    // create group specified by PATH
    hid_t group_id = H5Gcreate1(file_id, path, H5P_DEFAULT);
    H5Gclose(group_id);

    // initialize snapshot count
    hist->snapshot = 0;

    // store dimensionality
    hist->n_snapshots = n_snapshots;

    // save bins (include an extra rightmost bin??)
    hist->nbins1 = nbins1 + (adjust_right_bin ? 1 : 0);
    hist->nbins2 = nbins2 + (adjust_right_bin ? 1 : 0);

    // is this a 2-d histogram??
    hist->two_dim = (ylim != NULL);

    // create path to X
    GString *path_x = g_string_new(path);
    g_string_append_path(path_x, "x");
    // create path to Y
    GString *path_y = g_string_new(path);
    g_string_append_path(path_y, "y");


    if (hist->two_dim) {

	// === 2-D HISTOGRAM ===

	// create path to z
	GString *path_z = g_string_new(path);
	g_string_append_path(path_z, "z");

	// create the histogram
	hist->histogram = gsl_histogram2d_alloc(hist->nbins1, hist->nbins2);
	gsl_histogram2d_set_ranges_uniform(hist->histogram,
					   xlim[0], xlim[1] + (adjust_right_bin ? xlim[1] / nbins1 : 0),
					   ylim[0], ylim[1] + (adjust_right_bin ? ylim[1] / nbins2 : 0));

	// create, write and finalize the X dataset
	hdf_dataset_initialize(&hist->dataset_x, file_id, path_x->str, 1, &hist->nbins1, Q_HDF_TYPE_DOUBLE);
	hdf_dataset_write_all(&hist->dataset_x, ((gsl_histogram2d *)hist->histogram)->xrange);
	hdf_dataset_finalize(&hist->dataset_x);
	// create, write and finalize the Y dataset
	hdf_dataset_initialize(&hist->dataset_y, file_id, path_y->str, 1, &hist->nbins2, Q_HDF_TYPE_DOUBLE);
	hdf_dataset_write_all(&hist->dataset_y, ((gsl_histogram2d *)hist->histogram)->yrange);
	hdf_dataset_finalize(&hist->dataset_y);

	// create the Z (2-d) dataset
	hsize_t dims_z[] = {n_snapshots, hist->nbins1, hist->nbins2};
	hdf_dataset_initialize(&hist->dataset_z, file_id, path_z->str, 3, dims_z, Q_HDF_TYPE_DOUBLE);

	g_string_free(path_z, TRUE);

	
    } else {

	// === 1-D HISTOGRAM ===

	// create the histogram
	hist->histogram = gsl_histogram_alloc(hist->nbins1);
	gsl_histogram_set_ranges_uniform(hist->histogram, xlim[0], xlim[1] + (adjust_right_bin ? xlim[1] / nbins1 : 0));
	//hist->histogram = NULL;

	// create, write and finalize the X dataset
	hdf_dataset_initialize(&hist->dataset_x, file_id, path_x->str, 1, &hist->nbins1, Q_HDF_TYPE_DOUBLE);
	hdf_dataset_write_all(&hist->dataset_x, ((gsl_histogram *)hist->histogram)->range);
	hdf_dataset_finalize(&hist->dataset_x);

	// create the Y (2-d) dataset
	hsize_t dims_y[] = {n_snapshots, hist->nbins1};
	hdf_dataset_initialize(&hist->dataset_y, file_id, path_y->str, 2, dims_y, Q_HDF_TYPE_DOUBLE);
	
    }


    g_string_free(path_x, TRUE);
    g_string_free(path_y, TRUE);
    
}







// close/free all resources
void hdf_histogram_finalize(hdf_histogram_t *hist) {

    if (hist->snapshot < hist->n_snapshots)
	hdf_histogram_write(hist);

    // truncate and free the histograms (acc. to written snapshots)
    hsize_t tr_size[] = {hist->snapshot, hist->nbins1, hist->nbins2};

    if (hist->two_dim) {
	hdf_dataset_truncate(&hist->dataset_z, tr_size);
	hdf_dataset_finalize(&hist->dataset_z);
	gsl_histogram2d_free(hist->histogram);

    } else {
	hdf_dataset_truncate(&hist->dataset_y, tr_size);
	hdf_dataset_finalize(&hist->dataset_y);
	gsl_histogram_free(hist->histogram);
    }

    
}









// write the element and reset histogram
void hdf_histogram_write(hdf_histogram_t *hist) {

    hsize_t offsets[] = {hist->snapshot, 0, 0};
    hsize_t counts[] = {1, 0, 0};

    if (hist->two_dim) {

	// === 2-D HISTOGRAM ===
	counts[1] = ((gsl_histogram2d *)hist->histogram)->nx;
	counts[2] = ((gsl_histogram2d *)hist->histogram)->ny;
	hdf_dataset_write(&hist->dataset_z, ((gsl_histogram2d *)hist->histogram)->bin, offsets, counts);

	// reset the histogram
	gsl_histogram2d_reset(hist->histogram);

    } else {

	// === 1-D HISTOGRAM ===
	counts[1] = ((gsl_histogram *)hist->histogram)->n;
	hdf_dataset_write(&hist->dataset_y, ((gsl_histogram *)hist->histogram)->bin, offsets, counts);

	// reset the histogram
	gsl_histogram_reset(hist->histogram);
	
    }

    // advance snapshot
    ++ hist->snapshot;
    
}









// add the value to the histogram
void hdf_histogram_add(hdf_histogram_t *hist, double val) {

    gsl_histogram_increment(hist->histogram, val);
    
}




// add the value to the histogram
void hdf_histogram_add2(hdf_histogram_t *hist, double xval, double yval, double weight) {

    gsl_histogram2d_accumulate(hist->histogram, xval, yval, weight);
    
}






//======================================================================================
//======================================================================================




// create the objects

void hdf_histogram2_initialize(hdf_histogram2_t *hist, hid_t file_id, const char *path, 
			      int n_snapshots, int n_elements,
			      hsize_t nbins1, hsize_t nbins2,
			      double *xlim, double *ylim,
			      int adjust_right_bin)
{

    // create group specified by PATH
    hid_t group_id = H5Gcreate1(file_id, path, H5P_DEFAULT);
    H5Gclose(group_id);

    // initialize snapshot count
    hist->snapshot = 0;
    // initialize element within snapshot count
    hist->element = 0;

    // store dimensionality
    hist->n_snapshots = n_snapshots;
    hist->n_elements = n_elements;

    // save bins (include an extra rightmost bin??)
    hist->nbins1 = nbins1 + (adjust_right_bin ? 1 : 0);
    hist->nbins2 = nbins2 + (adjust_right_bin ? 1 : 0);

    

    // is this a 2-d histogram??
    hist->two_dim = (ylim != NULL);

    //printf("%s (nx=%d, ny=%d) %s\n", path, hist->nbins1, hist->nbins2, (hist->two_dim ? "2D" : "1D"));

    // create path to X
    GString *path_x = g_string_new(path);
    g_string_append_path(path_x, "x");
    // create path to Y
    GString *path_y = g_string_new(path);
    g_string_append_path(path_y, "y");


    if (hist->two_dim) {

	// === 2-D HISTOGRAM ===

	// create path to z
	GString *path_z = g_string_new(path);
	g_string_append_path(path_z, "z");

	// create the histogram
	hist->histogram = gsl_histogram2d_alloc(hist->nbins1, hist->nbins2);
	gsl_histogram2d_set_ranges_uniform(hist->histogram,
					   xlim[0], xlim[1] + (adjust_right_bin ? xlim[1] / nbins1 : 0),
					   ylim[0], ylim[1] + (adjust_right_bin ? ylim[1] / nbins2 : 0));

	// create, write and finalize the X dataset
	hdf_dataset_initialize(&hist->dataset_x, file_id, path_x->str, 1, &hist->nbins1, Q_HDF_TYPE_DOUBLE);
	hdf_dataset_write_all(&hist->dataset_x, ((gsl_histogram2d *)hist->histogram)->xrange);
	hdf_dataset_finalize(&hist->dataset_x);
	// create, write and finalize the Y dataset
	hdf_dataset_initialize(&hist->dataset_y, file_id, path_y->str, 1, &hist->nbins2, Q_HDF_TYPE_DOUBLE);
	hdf_dataset_write_all(&hist->dataset_y, ((gsl_histogram2d *)hist->histogram)->yrange);
	hdf_dataset_finalize(&hist->dataset_y);

	// create the Z (2-d) dataset
	hsize_t dims_z[] = {n_snapshots, n_elements, hist->nbins1, hist->nbins2};
	hdf_dataset_initialize(&hist->dataset_z, file_id, path_z->str, 4, dims_z, Q_HDF_TYPE_DOUBLE);

	g_string_free(path_z, TRUE);

	
    } else {

	// === 1-D HISTOGRAM ===

	// create the histogram
	hist->histogram = gsl_histogram_alloc(hist->nbins1);
	gsl_histogram_set_ranges_uniform(hist->histogram, xlim[0], xlim[1] + (adjust_right_bin ? xlim[1] / nbins1 : 0));
	//hist->histogram = NULL;

	// create, write and finalize the X dataset
	hdf_dataset_initialize(&hist->dataset_x, file_id, path_x->str, 1, &hist->nbins1, Q_HDF_TYPE_DOUBLE);
	hdf_dataset_write_all(&hist->dataset_x, ((gsl_histogram *)hist->histogram)->range);
	hdf_dataset_finalize(&hist->dataset_x);

	// create the Y (2-d) dataset
	hsize_t dims_y[] = {n_snapshots, n_elements, hist->nbins1};
	hdf_dataset_initialize(&hist->dataset_y, file_id, path_y->str, 3, dims_y, Q_HDF_TYPE_DOUBLE);
	
    }


    g_string_free(path_x, TRUE);
    g_string_free(path_y, TRUE);
    
}







// close/free all resources
void hdf_histogram2_finalize(hdf_histogram2_t *hist) {

    if (hist->snapshot < hist->n_snapshots)
	hdf_histogram2_write(hist);

    // truncate and free the histograms (acc. to written snapshots)
    hsize_t tr_size[] = {hist->snapshot, hist->n_elements, hist->nbins1, hist->nbins2};

    if (hist->two_dim) {
	hdf_dataset_truncate(&hist->dataset_z, tr_size);
	hdf_dataset_finalize(&hist->dataset_z);
	gsl_histogram2d_free(hist->histogram);

    } else {
	hdf_dataset_truncate(&hist->dataset_y, tr_size);
	hdf_dataset_finalize(&hist->dataset_y);
	gsl_histogram_free(hist->histogram);
    }

    
}









// write the element and reset histogram
void hdf_histogram2_write(hdf_histogram2_t *hist) {

    hsize_t offsets[] = {hist->snapshot, hist->element, 0, 0};
    hsize_t counts[] = {1, 1, 0, 0};

    if (hist->two_dim) {

	// === 2-D HISTOGRAM ===
	counts[2] = ((gsl_histogram2d *)hist->histogram)->nx;
	counts[3] = ((gsl_histogram2d *)hist->histogram)->ny;
	hdf_dataset_write(&hist->dataset_z, ((gsl_histogram2d *)hist->histogram)->bin, offsets, counts);

	// reset the histogram
	gsl_histogram2d_reset(hist->histogram);

    } else {

	// === 1-D HISTOGRAM ===
	counts[2] = ((gsl_histogram *)hist->histogram)->n;
	hdf_dataset_write(&hist->dataset_y, ((gsl_histogram *)hist->histogram)->bin, offsets, counts);

	// reset the histogram
	gsl_histogram_reset(hist->histogram);
	
    }

    // increment element
    ++ hist->element;

    // advance snapshot??
    if (hist->element == hist->n_elements) {
	++ hist->snapshot;
	hist->element = 0;
    }
    
}





// increment the hist's snapshot (and reset element) -- NO WRITE WILL BE PERFORMED
void hdf_histogram2_advance(hdf_histogram2_t *hist) {

    ++ hist->snapshot;
    hist->element = 0;
    
}









// add the value to the histogram
void hdf_histogram2_add(hdf_histogram2_t *hist, double val) {

    gsl_histogram_increment(hist->histogram, val);
    
}




// add the value to the histogram
void hdf_histogram2_add2(hdf_histogram2_t *hist, double xval, double yval, double weight) {

    gsl_histogram2d_accumulate(hist->histogram, xval, yval, weight);
    
}

