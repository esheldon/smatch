#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "stack.h"
#include "tree.h"

void tree_insert(struct tree_node ** self, int64_t val, size_t index) {
    if(!(*self)) {
        *self = malloc(sizeof(struct tree_node));
        (*self)->left = (*self)->right = NULL;
        (*self)->val=val;
        (*self)->indices = szstack_new(1);
        szstack_push((*self)->indices, index);
        return;
    }
    if(val < (*self)->val) {
        tree_insert(&(*self)->left, val, index);
    } else if(val > (*self)->val) {
        tree_insert(&(*self)->right, val, index);
    } else {
        szstack_push((*self)->indices, index);
    }
}

struct tree_node* tree_delete(struct tree_node* self) {
    if (self != NULL) {
        self->indices = szstack_delete(self->indices);
        self->left    = tree_delete(self->left);
        self->right   = tree_delete(self->right);
        free(self);
    }
    return NULL;
}

struct tree_node* tree_find(struct tree_node* self, int64_t val) {
    if (!self) {
        return NULL;
    }
    if(val < self->val) {
        return tree_find(self->left, val);
    } else if(val > self->val) {
        return tree_find(self->right, val);
    } else {
        return self;
    }
}

void tree_print_padding( char ch, int n )
{
    int i;

    for ( i = 0; i < n; i++ ) {
        fprintf(stderr,"%c", ch);
    }
}

// send level=0 at first
void tree_print( struct tree_node *self, int level )
{
    if ( self == NULL ) {
        tree_print_padding ( '\t', level );
        fprintf(stderr,"~");
    } else {
        tree_print( self->right, level + 1 );
        tree_print_padding ( '\t', level );
        fprintf(stderr,"%ld\n", self->val);
        tree_print( self->left, level + 1 );
    }
}





