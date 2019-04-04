/*
   catalog points carry extra information, such as
   radius and cos radius
*/
#ifndef _CATPOINT_H
#define _CATPOINT_H

#include <stdint.h>

typedef struct {
    double x;
    double y;
    double z;
    double radius; // radians
    double cos_radius;
} CatPoint;

#endif
