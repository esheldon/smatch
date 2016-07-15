#ifndef _POINT_H
#define _POINT_H

#include <stdint.h>

typedef struct {
    double x;
    double y;
    double z;
    double radius; // radians
    double cos_radius;
} Point;

#endif
