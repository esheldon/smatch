#ifndef _TREE_HEADER
#define _TREE_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "stack.h"

struct tree_node {
    int64_t val;
    struct szstack* indices;
    struct tree_node* right, * left;
};


void tree_insert(struct tree_node ** self, int64_t val, size_t index);
struct tree_node* tree_find(struct tree_node* self, int64_t val);


void tree_print( struct tree_node *self, int level );
void tree_print_padding( char ch, int n );

struct tree_node* tree_delete(struct tree_node* self);
#endif
