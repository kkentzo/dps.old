/*

===========
bitvector.h
===========

A struct and functions to manipulate a vector of bits of a specified
length. Functionality includes accessing and mutating (flipping)
individual bits within the vector.

A note regarding the string representation (description): the string
should be read from right to left, so the first bit is the rightmost
character and the last bit the leftmost character.

 */


// ==========================================================
//              UTILITY STRUCTS AND FUNCTIONS
// ==========================================================

typedef struct {
    int index;
    int offset;
} pos_t;


// size_of_length : returns the number of chars needed to store BITS bits
int size_of_length(int bits);

// return the pos_t struct (index , offset) for bit position BIT
void find_pos(pos_t *pos, int bit);

// returns the hamming distance between two chars
int hamming(char c1, char c2);

// counting
int count_ones(char c);



//==================================================================
typedef struct {

    char *buf; // an array of chars
    int buf_size; // the number of chars in the data array
    int len; // the number of bits stored in the buffer

} bitvector_t;

//==================================================================
// allocate the bitvector (all bits are set to zero)
bitvector_t *bv_new( int length);

// create a new bitvector from a string
bitvector_t *bv_new_from_string(const char *bv_string);

// return an exact (deep) copy of the bitvector
// don't forget to free it!
bitvector_t *bv_copy(bitvector_t *bv);

// free the bitvector
void bv_free(bitvector_t *bv);

// reset the bitvector to a random sequence of bits 
void bv_randomize(bitvector_t *bv, gsl_rng *rng);

// set all bits on
void bv_set_all_on(bitvector_t *bv);
// set all bits off
void bv_set_all_off(bitvector_t *bv);

// check whether the bit at a specified position is on (1) or off (0)
gboolean bv_is_on(bitvector_t *bv, int n);
gboolean bv_is_off(bitvector_t *bv, int n);

// return a count of 1's (or 0's) in the bitvector
int bv_count_ones(bitvector_t *bv);
int bv_count_zeros(bitvector_t *bv);

// turn on (or off) a bit at a specified position
void bv_set_on(bitvector_t *bv, int n);
void bv_set_off(bitvector_t *bv, int n);

// flip the bit at the specified position
void bv_toggle(bitvector_t *bv, int n);
// flip n bits at random
void bv_rtoggle(bitvector_t *bv, int n, gsl_rng *rng);

// flip all bits in the bitvector
void bv_complement(bitvector_t *bv);

// return the hamming distance between two bitvectors
int bv_hamming(bitvector_t *bv1, bitvector_t *bv2);

// returns true if the two bitvectors are equal
gboolean bv_equal(bitvector_t *bv1, bitvector_t *bv2);

// STRING REPRESENTATIONS

// return a string with the bitvector's string representation
GString *bv_description(bitvector_t *bv);

// append the bitvector's string representation to the supplied string
void bv_append_description(bitvector_t *bv, GString *st);

// write a string representation of the bitvector in st
// the GString is first emptied of anything it might hold
void bv_update_description(bitvector_t *bv, GString *st);

// print a representation of the bitstring on screen
void bv_print(bitvector_t *bv);
