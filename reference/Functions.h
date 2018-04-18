#ifndef FUNCTIONS_H
#define FUNCTIONS_H

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T> int nsgn(T val) {
    return (val < T(0) ? -1 : 1);
}

double length(double x, double y);
double cubicPulse(double x);
void rotate(double &x, double &y, double phi);
double triangleOccupancy(double out1, double in, double out2);
double trapezoidOccupancy(double out1, double out2, double in1, double in2);
double occupancy(double d11, double d12, double d21, double d22);

enum CellType {
    CELL_FLUID,
    CELL_SOLID
};

#endif
