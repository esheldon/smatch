#ifndef _CAT_H
#define _CAT_H

#include <stdint.h>
#include "vector.h"
#include "catpoint.h"

typedef struct {
    CatPoint point;
    // remove
    match_vector *matches;
    // list of healpix ids that intersect the disc around
    // this point
    lvector* disc_pixels;
} CatalogEntry;

typedef struct {
    size_t size;
    CatalogEntry* data;
} Catalog;

Catalog* cat_new(size_t n);
void cat_free_data(Catalog *self);

#define cat_free(cat) do {                                                   \
    if ((cat)) {                                                             \
        cat_free_data(cat);                                                  \
        free((cat)->data);                                                   \
        free((cat));                                                         \
        (cat)=NULL;                                                          \
    }                                                                        \
} while(0)

#endif
