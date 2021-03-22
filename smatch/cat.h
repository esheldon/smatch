#ifndef _CAT_H
#define _CAT_H

#include <stdint.h>
#include "vector.h"
#include "catpoint.h"

typedef struct {
    CatPoint point;

    // to hold matches
    match_vector *matches;

    // list of healpix ids that intersect the disc around
    // this point
    lvector* disc_pixels;

} CatalogEntry;

// create a catalog entry, including making
// the match and pixels vectors

CatalogEntry* cat_entry_new(void);

#define cat_entry_free(entry) do {                                           \
    if ((entry)) {                                                           \
        vector_free(entry->disc_pixels);                                     \
        vector_free(entry->matches);                                         \
        free((entry));                                                       \
        (entry)=NULL;                                                        \
    }                                                                        \
} while(0)


/*

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
*/

#endif
