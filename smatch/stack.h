// This header was auto-generated
#ifndef _STACK_H
#define _STACK_H
#include <stdint.h>

#define STACK_PUSH_REALLOC_MULT 1
#define STACK_PUSH_REALLOC_MULTVAL 1.5
#define STACK_PUSH_INITSIZE 50

#ifndef float64
#define float64 double
#endif

struct i64stack {
    size_t size;            // number of elements that are visible to the user
    size_t allocated_size;  // number of allocated elements in data vector
    size_t push_realloc_style; // Currently always STACK_PUSH_REALLOC_MULT, 
                               // which is reallocate to allocated_size*realloc_multval
    size_t push_initsize;      // default size on first push, default STACK_PUSH_INITSIZE 
    double realloc_multval; // when allocated size is exceeded while pushing, 
                            // reallocate to allocated_size*realloc_multval, default 
                            // STACK_PUSH_REALLOC_MULTVAL
                            // if allocated_size was zero, we allocate to push_initsize
    int64_t* data;
};

struct i64stack* i64stack_new(size_t num);

// if size > allocated size, then a reallocation occurs
// if size <= internal size, then only the ->size field is reset
// use i64stack_realloc() to reallocate the data vector and set the ->size
void i64stack_resize(struct i64stack* stack, size_t newsize);

// perform reallocation on the underlying data vector. This does
// not change the size field unless the new size is smaller
// than the viewed size
void i64stack_realloc(struct i64stack* stack, size_t newsize);

// completely clears memory in the data vector
void i64stack_clear(struct i64stack* stack);

// clears all memory and sets pointer to NULL
// usage: stack=i64tack_delete(stack);
struct i64stack* i64stack_delete(struct i64stack* stack);

// if reallocation is needed, size is increased by 50 percent
// unless size is zero, when it 100 are allocated
void i64stack_push(struct i64stack* stack, int64_t val);
// pop the last element and decrement size; no reallocation is performed
// if empty, INT64_MIN is returned
int64_t i64stack_pop(struct i64stack* stack);

int __i64stack_compare_el(const void *a, const void *b);
void i64stack_sort(struct i64stack* stack);
int64_t* i64stack_find(struct i64stack* stack, int64_t el);

struct f64stack {
    size_t size;            // number of elements that are visible to the user
    size_t allocated_size;  // number of allocated elements in data vector
    size_t push_realloc_style; // Currently always STACK_PUSH_REALLOC_MULT, 
                               // which is reallocate to allocated_size*realloc_multval
    size_t push_initsize;      // default size on first push, default STACK_PUSH_INITSIZE 
    double realloc_multval; // when allocated size is exceeded while pushing, 
                            // reallocate to allocated_size*realloc_multval, default 
                            // STACK_PUSH_REALLOC_MULTVAL
                            // if allocated_size was zero, we allocate to push_initsize
    float64* data;
};

struct f64stack* f64stack_new(size_t num);

// if size > allocated size, then a reallocation occurs
// if size <= internal size, then only the ->size field is reset
// use f64stack_realloc() to reallocate the data vector and set the ->size
void f64stack_resize(struct f64stack* stack, size_t newsize);

// perform reallocation on the underlying data vector. This does
// not change the size field unless the new size is smaller
// than the viewed size
void f64stack_realloc(struct f64stack* stack, size_t newsize);

// completely clears memory in the data vector
void f64stack_clear(struct f64stack* stack);

// clears all memory and sets pointer to NULL
// usage: stack=f64tack_delete(stack);
struct f64stack* f64stack_delete(struct f64stack* stack);

// if reallocation is needed, size is increased by 50 percent
// unless size is zero, when it 100 are allocated
void f64stack_push(struct f64stack* stack, float64 val);
// pop the last element and decrement size; no reallocation is performed
// if empty, INT64_MIN is returned
float64 f64stack_pop(struct f64stack* stack);

int __f64stack_compare_el(const void *a, const void *b);
void f64stack_sort(struct f64stack* stack);
float64* f64stack_find(struct f64stack* stack, float64 el);

struct szstack {
    size_t size;            // number of elements that are visible to the user
    size_t allocated_size;  // number of allocated elements in data vector
    size_t push_realloc_style; // Currently always STACK_PUSH_REALLOC_MULT, 
                               // which is reallocate to allocated_size*realloc_multval
    size_t push_initsize;      // default size on first push, default STACK_PUSH_INITSIZE 
    double realloc_multval; // when allocated size is exceeded while pushing, 
                            // reallocate to allocated_size*realloc_multval, default 
                            // STACK_PUSH_REALLOC_MULTVAL
                            // if allocated_size was zero, we allocate to push_initsize
    size_t* data;
};

struct szstack* szstack_new(size_t num);

// if size > allocated size, then a reallocation occurs
// if size <= internal size, then only the ->size field is reset
// use szstack_realloc() to reallocate the data vector and set the ->size
void szstack_resize(struct szstack* self, size_t newsize);

// perform reallocation on the underlying data vector. This does
// not change the size field unless the new size is smaller
// than the viewed size
void szstack_realloc(struct szstack* self, size_t newsize);

// completely clears memory in the data vector
void szstack_clear(struct szstack* self);

// clears all memory and sets pointer to NULL
// usage: stack=sztack_delete(stack);
struct szstack* szstack_delete(struct szstack* self);

// if reallocation is needed, size is increased by 50 percent
// unless size is zero, when it 100 are allocated
void szstack_push(struct szstack* self, size_t val);
// pop the last element and decrement size; no reallocation is performed
// if empty, INT64_MIN is returned
size_t szstack_pop(struct szstack* self);

int __szstack_compare_el(const void *a, const void *b);
void szstack_sort(struct szstack* self);
size_t* szstack_find(struct szstack* self, size_t el);

#endif  // header guard
