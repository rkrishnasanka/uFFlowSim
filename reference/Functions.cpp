//#include "Functions.h"

#define _USE_MATH_DEFINES

#include <algorithm>
#include <stdint.h>
#include <math.h>

using namespace std;

double length(double x, double y) {
    return sqrt(x*x + y*y);
}

double cubicPulse(double x) {
    x = min(fabs(x), 1.0);
    return 1.0 - x*x*(3.0 - 2.0*x);
}

void rotate(double &x, double &y, double phi) {
    double tmpX = x, tmpY = y;
    x =  cos(phi)*tmpX + sin(phi)*tmpY;
    y = -sin(phi)*tmpX + cos(phi)*tmpY;
}

double triangleOccupancy(double out1, double in, double out2) {
	return 0.5*in*in/((out1 - in)*(out2 - in));
}

double trapezoidOccupancy(double out1, double out2, double in1, double in2) {
	return 0.5*(-in1/(out1 - in1) - in2/(out2 - in2));
}

double occupancy(double d11, double d12, double d21, double d22) {
	double ds[] = {d11, d12, d22, d21};

    /* Compute mask */
	uint8_t b = 0;
	for (int i = 3; i >= 0; i--)
		b = (b << 1) | (ds[i] < 0.0 ? 1 : 0);
    
	switch (b) {
    /* All outside */
	case 0x0: return 0.0;
    /* One inside */
	case 0x1: return triangleOccupancy(d21, d11, d12);
	case 0x2: return triangleOccupancy(d11, d12, d22);
	case 0x4: return triangleOccupancy(d12, d22, d21);
	case 0x8: return triangleOccupancy(d22, d21, d11);
    /* One outside */
	case 0xE: return 1.0 - triangleOccupancy(-d21, -d11, -d12);
	case 0xD: return 1.0 - triangleOccupancy(-d11, -d12, -d22);
	case 0xB: return 1.0 - triangleOccupancy(-d12, -d22, -d21);
	case 0x7: return 1.0 - triangleOccupancy(-d22, -d21, -d11);
    /* Two adjacent inside */
	case 0x3: return trapezoidOccupancy(d21, d22, d11, d12);
	case 0x6: return trapezoidOccupancy(d11, d21, d12, d22);
	case 0x9: return trapezoidOccupancy(d12, d22, d11, d21);
	case 0xC: return trapezoidOccupancy(d11, d12, d21, d22);
    /* Two opposed inside */
	case 0x5: return triangleOccupancy(d11, d12, d22) +
		             triangleOccupancy(d22, d21, d11);
	case 0xA: return triangleOccupancy(d21, d11, d12) +
		             triangleOccupancy(d12, d22, d21);
    /* All inside */
	case 0xF: return 1.0;
	}
    
    return 0.0;
}
