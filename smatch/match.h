#ifndef _MATCH_H
#define _MATCH_H

#include <stdint.h>

#define MATCHSTACK_PUSH_REALLOC_MULT 1
#define MATCHSTACK_PUSH_REALLOC_MULTVAL 1.5
#define MATCHSTACK_PUSH_INITSIZE 10

struct match {
    size_t cat_ind;
    size_t input_ind;
    double cosdist;
};

struct matchstack {

    size_t size;            // number of elements that are visible to the user
    size_t allocated_size;  // number of allocated elements in data vector
    size_t push_realloc_style; // Currently always STACK_PUSH_REALLOC_MULT, 
                               // which is reallocate to allocated_size*realloc_multval
    size_t push_initsize;      // default size on first push, default STACK_PUSH_INITSIZE 
    double realloc_multval; // when allocated size is exceeded while pushing, 
                            // reallocate to allocated_size*realloc_multval, default 
                            // STACK_PUSH_REALLOC_MULTVAL
                            // if allocated_size was zero, we allocate to push_initsize
    struct match* data;

};

struct matchstack* matchstack_new(void);
void matchstack_realloc(struct matchstack* ms, size_t newsize);
void matchstack_resize(struct matchstack* ms, size_t newsize);
void matchstack_push(struct matchstack* ms, 
                     size_t cat_ind, size_t input_ind, double cosdist);

int match_compare(const void *a, const void *b);
void matchstack_sort(struct matchstack* ms);
struct matchstack* matchstack_delete(struct matchstack* ms);
#endif
