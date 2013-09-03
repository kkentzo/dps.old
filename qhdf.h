/*

  qhdf.h

  Supplies classes and function to simplify working with :

  === 1/ HDF tables (see hdf_table_* functions) ===

      represents a table of records, each having N individual fields ,
      stored in memory (and on disk) as an HDF dataset in a file
      can be saved to and read from an HDF file

      ALL FIELDS ARE OF THE SAME TYPE
      this way, explicit structs are not needed -- just pass a vector with the field values


  === 2/ HDF datasets (see hdf_dataset_* functions) ===

      represents an HDF dataset of a fixed dimensionality


  === 3/ HDF Histograms (see hdf_histogram_* functions) ===

      represents a 1-D or 2-D histogram across many _snapshots_ (e.g. time points)

      *** also available is a histogram (see hdf_histogram2_t) that
      supports multiple elements (i.e. 1- or 2-d histograms) within a
      single snapshot. Imagine a 2-d plane where we're supposed to be
      stacking rectangular objects (the histograms) along the x-axis
      are the snapshots and along the y-axis are the elements.  we
      stack elements until the max height is reached and then we go to
      the next position (snapshot) and start stacking again. The
      hdf_histogram2_write() operation does exactly that. Use
      N_SNAPSHOTS and N_ELEMENTS to customize this behaviour

 */




//======================================================================================
//======================================================================================
#define TABLE_DEFAULT_CHUNK_SIZE 1000



typedef enum {

    Q_HDF_TYPE_INT,
    Q_HDF_TYPE_DOUBLE,

} Q_HDF_TYPE;




typedef struct {

    // the HDF5 file in which the table will be stored
    hid_t file_id;
    // the name of the dataset in which the table will be stored (will be created)
    const char *dataset_name; 
    // the number of fields per record (equal-sized fields)
    size_t nfields; 
    // the total size of a record (calculated)
    size_t record_size;

    // attributes for creating/appending to the table
    hid_t *field_types; // an array with the field types
    size_t *field_sizes; // an array with the field sizes
    size_t *field_offsets; // an array with the field offsets
    const char **field_names; // the field names

    GPtrArray *generated_names; // automatically generated field names 

    // chunk size 
    hsize_t _chunk_size_;

    // buffer of records
    void *buf;
    int buf_idx;
    // how many records to buffer before writing to the file 
    size_t buf_size; // default value: TABLE_DEFAULT_CHUNK_SIZE

} hdf_table_t;




// initialize a table (hdf_table_t struct already statically exists)
// FIELD_NAMES can be NULL, in which case names are auto-generated as letters of
// the alphabet; the case were nfields >= LETTERS_SIZE is also handled (numbers are
// appended to letters accordingly)
void hdf_table_initialize(hdf_table_t *table, hid_t file_id, const char *dataset_name, 
			  size_t nfields, Q_HDF_TYPE table_type, const char **field_names);

// finalize the table (free internal resources)
void hdf_table_finalize(hdf_table_t *table);

// append the supplied record to the table
// DATA should be of length table->nfields
void hdf_table_append_record(hdf_table_t *table, void *data);

// write the internal buffer to the table and reset the buffer index
void hdf_table_flush(hdf_table_t *table);



// THIS WORKS BUT NEEDS DEVELOPMENT -- SEE ACTUAL FUNCTION IN table.c
// ===>
// attach (write) an attribute to the table
//void table_attach_attribute(hdf_table_t *table);


// ===========================================================================
//        READ A TABLE FROM A FILE AND RETURN IT IN MATRIX FORM
// ===========================================================================
// loads the requested table from open file FILE_ID
// writes the number of records in NRECORDS
// returns a newly allocated buffer with the table's contents
// remember to FREE the buffer
void *hdf_table_read(hid_t file_id, const char *table_name, hsize_t *nrecords);









//======================================================================================
//======================================================================================


typedef struct {

    // the HDF5 file in which the table will be stored
    hid_t file_id;
    // the dataset
    hid_t dataset_id;
    // the data type
    hid_t dataset_type;

    hsize_t ndim;
    hsize_t *dims;

} hdf_dataset_t;





void hdf_dataset_initialize(hdf_dataset_t *dataset, hid_t file_id, const char *name,
			    hsize_t ndim, hsize_t *dims, Q_HDF_TYPE dataset_type);

void hdf_dataset_finalize(hdf_dataset_t *dataset);

void hdf_dataset_truncate(hdf_dataset_t *dataset, hsize_t *tr_size);

void hdf_dataset_write(hdf_dataset_t *dataset, void *data, hsize_t *offsets, hsize_t *counts);

void hdf_dataset_write_all(hdf_dataset_t *dataset, void *data);




//======================================================================================
//======================================================================================

// A HISTOGRAM
typedef struct {

    hdf_dataset_t dataset_x;
    hdf_dataset_t dataset_y;
    hdf_dataset_t dataset_z;

    int two_dim;

    int snapshot;

    int n_snapshots;

    hsize_t nbins1;
    hsize_t nbins2;

    void *histogram;

} hdf_histogram_t;


// create the histogram -- if YLIM is NULL then it's a 1-d histogram, otherwise 2-D
// an hdf group will be created acc to the string specified by PATH
// if ADJUST_RIGHT_BIN the histograms include an extra rightmost bin
// (it's all handled automatically - do not mess with it from the outside)
void hdf_histogram_initialize(hdf_histogram_t *hist, hid_t file_id, const char *path, 
			      int n_snapshots, hsize_t nbins1, hsize_t nbins2,
			      double *xlim, double *ylim,
			      int adjust_right_bin);

// close/free all resources
void hdf_histogram_finalize(hdf_histogram_t *hist);


// write the current SNAPSHOT and reset histogram to zeros
void hdf_histogram_write(hdf_histogram_t *hist);

// add the value to the histogram
void hdf_histogram_add(hdf_histogram_t *hist, double val);
void hdf_histogram_add2(hdf_histogram_t *hist, double xval, double yval, double weight);





//======================================================================================
//======================================================================================

// A HISTOGRAM WITH MULTIPLE ELEMENTS PER SNAPSHOT
typedef struct {

    hdf_dataset_t dataset_x;
    hdf_dataset_t dataset_y;
    hdf_dataset_t dataset_z;

    int two_dim;

    int snapshot;
    int element;

    int n_snapshots;
    int n_elements;

    hsize_t nbins1;
    hsize_t nbins2;

    void *histogram;

} hdf_histogram2_t;


// create the histogram -- if YLIM is NULL then it's a 1-d histogram, otherwise 2-D
// an hdf group will be created acc to the string specified by PATH
// if ADJUST_RIGHT_BIN the histograms include an extra rightmost bin
// (it's all handled automatically - do not mess with it from the outside)
void hdf_histogram2_initialize(hdf_histogram2_t *hist, hid_t file_id, const char *path, 
			       int n_snapshots, int n_elements,
			       hsize_t nbins1, hsize_t nbins2,
			       double *xlim, double *ylim,
			       int adjust_right_bin);

// close/free all resources
void hdf_histogram2_finalize(hdf_histogram2_t *hist);


// write the current ELEMENT and reset histogram to zeros
// advances the snapshot as well, if the written element is the max one
void hdf_histogram2_write(hdf_histogram2_t *hist);

// increment the hist's snapshot (and reset element) -- NO WRITE WILL BE PERFORMED
void hdf_histogram2_advance(hdf_histogram2_t *hist);


// add the value to the histogram
void hdf_histogram2_add(hdf_histogram2_t *hist, double val);
void hdf_histogram2_add2(hdf_histogram2_t *hist, double xval, double yval, double weight);
