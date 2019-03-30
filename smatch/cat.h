#ifndef _CAT_H
#define _CAT_H

#include <stdint.h>
#include "vector.h"
#include "point.h"

typedef struct {
    Point point;
    match_vector *matches;
} CatalogEntry;

typedef struct {
    size_t size;
    CatalogEntry* data;
} Catalog;

Catalog* cat_new(size_t n);
void cat_free_matches(Catalog *self);

#define cat_free(cat) do {                                                   \
    if ((cat)) {                                                             \
        cat_free_matches(cat);                                               \
        free((cat)->data);                                                   \
        free((cat));                                                         \
        (cat)=NULL;                                                          \
    }                                                                        \
} while(0)

#endif
