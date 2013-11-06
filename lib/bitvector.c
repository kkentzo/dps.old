#include <math.h>
#include <assert.h>
#include <string.h>
#include <glib.h>
#include <gsl/gsl_rng.h>

#include "bitvector.h"
#include "sampling.h"





// ==========================================================
//                 UTILITY FUNCTIONS
// ==========================================================

// returns the number of char bytes needed to store BITS bits
int size_of_length(int bits){

    int size = bits / CHAR_BIT;

    if (bits % CHAR_BIT)
	size += 1;
    
    return size;
		    
}



// ==========================================================
// returns the index and offset of the specified bit
// inside the char vector
void find_pos(pos_t *pos, int n) {

    div_t result = div(n, CHAR_BIT);

    pos->index = result.quot;
    pos->offset = result.rem;

}



// ==========================================================
// returns the haming distance between the two chars
int hamming(char c1, char c2) {
    
    char xor_val = c1 ^ c2;
    int dist = 0;

    while (xor_val) {
	dist++;
	xor_val &= xor_val - 1;
    }

    return dist;
	    
}


// ==========================================================
// counts the number of ones in the char
int count_ones(char c) {

    int i, cnt = 0;

    for (i=0; i<CHAR_BIT; i++)
	if (c & 1 << i)
	    ++cnt;


    return cnt;

}




// ==========================================================
//                 BITVECTOR FUNCTIONS
// ==========================================================
bitvector_t *bv_new(int length) {

    assert(length > 0);

    bitvector_t *bv = malloc(sizeof(bitvector_t));
    // store number of bits
    bv->len = length;
    // calculate buffer size
    bv->buf_size = size_of_length(length);
    // allocate space for buffer (it's zeroed)
    bv->buf = (char *) calloc(bv->buf_size, sizeof(char));

    return bv;
    
}





// create a new bitvector from a string
bitvector_t *bv_new_from_string(const char *bv_string) {

    // create a gstring
    GString *st = g_string_new(bv_string);
    // create the bitvector
    bitvector_t *bv = bv_new(st->len);
    // read the string into the bitvector
    int i;
    for (i=0; i<st->len; i++)
	if (bv_string[i] == '1')
	    bv_set_on(bv, i);
    
    g_string_free(st, TRUE);

    return bv;
    
}





bitvector_t *bv_copy(bitvector_t *bv) {

    // allocate space for copy
    bitvector_t *bv1 = bv_new(bv->len);
    // copy the buffer
    memmove(bv1->buf, bv->buf, bv->buf_size);
    return bv1;
    
}



void bv_free(bitvector_t *bv) {

    // free the buffer
    free(bv->buf);
    // free the struct
    free(bv);

}



// reset the bitvector to a random sequence of bits 
void bv_randomize(bitvector_t *bv, gsl_rng *rng) {

    int i;
    for (i=0; i<bv->len; i++)
	if (gsl_rng_uniform(rng) > 0.5)
	    bv_set_on(bv, i);
	else
	    bv_set_off(bv, i);
    
}





// set all bits on
void bv_set_all_on(bitvector_t *bv) {

    // fill bitvector with ones
    memset(bv->buf, ~0, bv->buf_size);
    
}




// set all bits off
void bv_set_all_off(bitvector_t *bv) {

    // fill bitvector with zeros
    memset(bv->buf, 0, bv->buf_size);
    
}





gboolean bv_is_on(bitvector_t *bv, int n) {

    // get position of n^th bit
    pos_t pos;
    find_pos(&pos, n);
    return (*(bv->buf + pos.index) & 1 << pos.offset) > 0;
}



gboolean bv_is_off(bitvector_t *bv, int n) {

    return ! bv_is_on(bv, n);
    
}



int bv_count_ones(bitvector_t *bv) {

    int i, cnt = 0;
    for (i=0; i<bv->buf_size; i++)
	cnt += count_ones(*(bv->buf + i));
    return cnt;

}


int bv_count_zeros(bitvector_t *bv) {

    return bv->len - bv_count_ones(bv);
    
}



void bv_set_on(bitvector_t *bv, int n) {

    pos_t pos;
    if (n < bv->len) {
	find_pos(&pos, n);
	*(bv->buf + pos.index) |= 1 << pos.offset;
    }
    
}


void bv_set_off(bitvector_t *bv, int n) {

    pos_t pos;
    if (n < bv->len) {
	find_pos(&pos, n);
	*(bv->buf + pos.index) &= ~ (1 << pos.offset);
    }
    
}




void bv_complement(bitvector_t *bv) {

    int i;
    for (i=0; i<bv->len; i++)
	bv_toggle(bv, i);
    
}


void bv_toggle(bitvector_t *bv, int n) {

    pos_t pos;
    if (n < bv->len) {
	find_pos(&pos, n);
	*(bv->buf + pos.index) ^= 1 << pos.offset;
    }
    
}



void bv_rtoggle(bitvector_t *bv, int n, gsl_rng *rng) {

    int i, positions[n];
    // select n random unique positions in the bitvector
    // ranging in [0, bv->len)
    sample(n, bv->len, positions, NULL, FALSE, FALSE, rng);
    //uselect(positions, n, 0, bv->len, rng);
    // flip selected bits
    for (i=0; i<n; i++)
	bv_toggle(bv, positions[i]);
    
}



int bv_hamming(bitvector_t *bv1, bitvector_t *bv2) {

    // assert the two bitvectors have the same length
    assert( bv1->len == bv2->len );

    int i, distance=0;
    for (i=0; i<bv1->buf_size; i++)
	distance += hamming(*(bv1->buf + i), *(bv2->buf + i));

    return distance;
    
}



gboolean bv_equal(bitvector_t *bv1, bitvector_t *bv2) {

    return bv_hamming(bv1, bv2) == 0;
    
}



GString *bv_description(bitvector_t *bv) {

    GString *st = g_string_sized_new(bv->len);
    bv_append_description(bv, st);
    return st;
    
}


void bv_append_description(bitvector_t *bv, GString *st) {

    int i;
    for (i=0; i<bv->len; i++)
	g_string_append_c(st, (bv_is_on(bv, i) ? '1' : '0'));
    
}


void bv_update_description(bitvector_t *bv, GString *st) {

    // empty string
    g_string_assign(st, "");
    // append bitvector to string
    bv_append_description(bv, st);

}



void bv_print(bitvector_t *bv) {

    GString *st = g_string_new("");
    bv_append_description(bv, st);
    printf("%s", st->str);
    g_string_free(st, TRUE);
    
}
