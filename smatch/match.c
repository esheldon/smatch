#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "match.h"

struct match* match_new(size_t cat_ind, size_t input_ind, double cosdist) {
    struct match* match=NULL;

    match=malloc(sizeof(struct match));
    if (match==NULL) {
        fprintf(stderr,"Error allocating match\n");
        return NULL;
    }

    match->cat_ind=cat_ind;
    match->input_ind=input_ind;
    match->cosdist=cosdist;

    return match;
}

struct matchstack* matchstack_new(void) {
    struct matchstack* ms=NULL;

    ms=malloc(sizeof(struct matchstack));
    if (ms==NULL) {
        fprintf(stderr,"Error allocating match\n");
        return NULL;
    }

    ms->size=0;
    ms->data=NULL;

    ms->push_realloc_style = MATCHSTACK_PUSH_REALLOC_MULT;
    ms->push_initsize = MATCHSTACK_PUSH_INITSIZE;
    ms->realloc_multval = MATCHSTACK_PUSH_REALLOC_MULTVAL;

    return ms;
}


void matchstack_realloc(struct matchstack* ms, size_t newsize) {

    size_t oldsize = ms->allocated_size;
    if (newsize != oldsize) {
        struct match* newdata=NULL;
        size_t elsize = sizeof(struct match);

        newdata = realloc(ms->data, newsize*elsize);
        if (newdata == NULL) {
            fprintf(stderr,"failed to reallocate\n");
            return;
        }

        if (newsize > ms->allocated_size) {
            // the allocated size is larger.  make sure to initialize the new
            // memory region.  This is the area starting from index [oldsize]
            size_t num_new_bytes = (newsize-oldsize)*elsize;
            memset(&newdata[oldsize], 0, num_new_bytes);
        } else if (ms->size > newsize) {
            // The viewed size is larger than the allocated size in this case,
            // we must set the size to the maximum it can be, which is the
            // allocated size
            ms->size = newsize;
        }

        ms->data = newdata;
        ms->allocated_size = newsize;
    }

}

void matchstack_resize(struct matchstack* ms, size_t newsize) {
   if (newsize > ms->allocated_size) {
       matchstack_realloc(ms, newsize);
   }

   ms->size = newsize;
}

void matchstack_push(struct matchstack* ms, size_t cat_ind, size_t input_ind, double cosdist) {
    // see if we have already filled the available data vector
    // if so, reallocate to larger storage
    struct match* m;
    if (ms->size == ms->allocated_size) {

        size_t newsize;
        if (ms->allocated_size == 0) {
            newsize=ms->push_initsize;
        } else {
            // currenly we always use the multiplier reallocation  method.
            if (ms->push_realloc_style != MATCHSTACK_PUSH_REALLOC_MULT) {
                fprintf(stderr,"Currently only support push realloc style MATCHSTACK_PUSH_REALLOC_MULT\n");
                exit(EXIT_FAILURE);
            }
            // this will "floor" the size
            newsize = (size_t)(ms->allocated_size*ms->realloc_multval);
            // we want ceiling
            newsize++;
        }

        matchstack_realloc(ms, newsize);

    }

    ms->size++;

    m = &ms->data[ms->size-1];
    m->cat_ind = cat_ind;
    m->input_ind = input_ind;
    m->cosdist = cosdist;

}

struct matchstack* matchstack_delete(struct matchstack* ms) {
    if (ms != NULL) {
        if (ms->data != NULL) {
            free(ms->data);
        }
        free(ms);
    }

    return NULL;
}
