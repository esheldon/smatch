#include <stdlib.h>
#include <stdio.h>
#include "vector.h"
#include "cat.h"

Catalog* cat_new(size_t n)
{
    size_t i=0;
    Catalog* self = calloc(1,sizeof(Catalog));

    if (self == NULL) {
        fprintf(stderr,"Could not allocate Catalog\n");
        return NULL;
    }

    self->size = n;
    self->data = calloc(self->size, sizeof(CatalogEntry));
    if (self->data == NULL) {
        fprintf(stderr,"Could not allocate data for Catalog\n");
        exit(1);
    }

    for (i=0; i<self->size; i++) {
        self->data[i].matches = match_vector_new();
        self->data[i].disc_pixels = lvector_new();
    }

    return self;
}

void cat_free_data(Catalog* self)
{
    size_t i=0;

    if (self) {
        for (i=0; i < self->size; i++) {
            //vector_free(self->data[i].matches);
            vector_free(self->data[i].disc_pixels);
        }
    }
}


