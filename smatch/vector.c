// This file was auto-generated using vectorgen
// most array methods are generic, see vector.h

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <float.h>
#include "vector.h"

point_vector* point_vector_new(void) {
    point_vector* self = calloc(1,sizeof(point_vector));
    if (self == NULL) {
        fprintf(stderr,"Could not allocate point_vector\n");
        return NULL;
    }

    self->capacity = VECTOR_INITCAP;

    self->data = calloc(self->capacity, sizeof(Point));
    if (self->data == NULL) {
        fprintf(stderr,"Could not allocate data for vector\n");
        exit(1);
    }

    return self;
}


point_vector* point_vector_copy(point_vector* self) {
    point_vector* vcopy=point_vector_new();
    vector_resize(vcopy, self->size);

    if (self->size > 0) {
        memcpy(vcopy->data, self->data, self->size*sizeof(Point));
    }

    return vcopy;
}

point_vector* point_vector_fromarray(Point* data, size_t size) {
    point_vector* self=point_vector_new();
    vector_resize(self, size);

    if (self->size > 0) {
        memcpy(self->data, data, size*sizeof(Point));
    }

    return self;
}

point_vector* point_vector_zeros(size_t num) {

    point_vector* self=point_vector_new();
    vector_resize(self, num);
    return self;
}


match_vector* match_vector_new(void) {
    match_vector* self = calloc(1,sizeof(match_vector));
    if (self == NULL) {
        fprintf(stderr,"Could not allocate match_vector\n");
        return NULL;
    }

    self->capacity = VECTOR_INITCAP;

    self->data = calloc(self->capacity, sizeof(Match));
    if (self->data == NULL) {
        fprintf(stderr,"Could not allocate data for vector\n");
        exit(1);
    }

    return self;
}


match_vector* match_vector_copy(match_vector* self) {
    match_vector* vcopy=match_vector_new();
    vector_resize(vcopy, self->size);

    if (self->size > 0) {
        memcpy(vcopy->data, self->data, self->size*sizeof(Match));
    }

    return vcopy;
}

match_vector* match_vector_fromarray(Match* data, size_t size) {
    match_vector* self=match_vector_new();
    vector_resize(self, size);

    if (self->size > 0) {
        memcpy(self->data, data, size*sizeof(Match));
    }

    return self;
}

match_vector* match_vector_zeros(size_t num) {

    match_vector* self=match_vector_new();
    vector_resize(self, num);
    return self;
}


lvector* lvector_new(void) {
    lvector* self = calloc(1,sizeof(lvector));
    if (self == NULL) {
        fprintf(stderr,"Could not allocate lvector\n");
        return NULL;
    }

    self->capacity = VECTOR_INITCAP;

    self->data = calloc(self->capacity, sizeof(int64_t));
    if (self->data == NULL) {
        fprintf(stderr,"Could not allocate data for vector\n");
        exit(1);
    }

    return self;
}


lvector* lvector_copy(lvector* self) {
    lvector* vcopy=lvector_new();
    vector_resize(vcopy, self->size);

    if (self->size > 0) {
        memcpy(vcopy->data, self->data, self->size*sizeof(int64_t));
    }

    return vcopy;
}

lvector* lvector_fromarray(int64_t* data, size_t size) {
    lvector* self=lvector_new();
    vector_resize(self, size);

    if (self->size > 0) {
        memcpy(self->data, data, size*sizeof(int64_t));
    }

    return self;
}

lvector* lvector_zeros(size_t num) {

    lvector* self=lvector_new();
    vector_resize(self, num);
    return self;
}


lvector* lvector_ones(size_t num) {

    lvector* self=lvector_new();
    for (size_t i=0; i<num; i++) {
        vector_push(self,1);
    }
    return self;
}

lvector* lvector_range(long min, long max) {

    lvector* self=lvector_new();
    for (long i=min; i<max; i++) {
        vector_push(self,i);
    }
    
    return self;
}

int __lvector_compare_el(const void *a, const void *b) {
    int64_t temp = 
        (  (int64_t) *( (int64_t*)a ) ) 
         -
        (  (int64_t) *( (int64_t*)b ) );
    if (temp > 0)
        return 1;
    else if (temp < 0)
        return -1;
    else
        return 0;
}


void lvector_sort(lvector* self) {
    qsort(self->data, self->size, sizeof(int64_t), __lvector_compare_el);
}
int64_t* lvector_find(lvector* self, int64_t el) {
    return (int64_t*) bsearch(&el, self->data, self->size, sizeof(int64_t), __lvector_compare_el);
}
