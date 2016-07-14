#ifndef _MATCH_H
#define _MATCH_H

#include <stdint.h>

typedef struct {
    size_t cat_ind;
    size_t input_ind;
    double cosdist;
} Match;


#endif
