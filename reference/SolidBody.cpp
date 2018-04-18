#define _USE_MATH_DEFINES

#include <algorithm>
#include <stdint.h>
#include <math.h>

#include "Functions.h"
#include "SolidBody.h"

using namespace std;

void SolidBody::globalToLocal(double &x, double &y) const {
    x -= _posX;
    y -= _posY;
    rotate(x, y, -_theta);
    x /= _scaleX;
    y /= _scaleY;
}

void SolidBody::localToGlobal(double &x, double &y) const {
    x *= _scaleX;
    y *= _scaleY;
    rotate(x, y, _theta);
    x += _posX;
    y += _posY;
}

double SolidBody::velocityX(double x, double y) const {
    return (_posY - y)*_velTheta + _velX;
}

double SolidBody::velocityY(double x, double y) const {
    return (x - _posX)*_velTheta + _velY;
}

void SolidBody::velocity(double &vx, double &vy, double x, double y) const {
    vx = velocityX(x, y);
    vy = velocityY(x, y);
}

void SolidBody::update(double timestep) {
    _posX  += _velX*timestep;
    _posY  += _velY*timestep;
    _theta += _velTheta*timestep;
}

double SolidBox::distance(double x, double y) const {
	x -= _posX;
	y -= _posY;
    rotate(x, y, -_theta);
	double dx = fabs(x) - _scaleX*0.5;
	double dy = fabs(y) - _scaleY*0.5;

	if (dx >= 0.0 || dy >= 0.0)
		return length(max(dx, 0.0), max(dy, 0.0));
	else
		return max(dx, dy);
}

void SolidBox::closestSurfacePoint(double &x, double &y) const {
	x -= _posX;
	y -= _posY;
	rotate(x, y, -_theta);
	double dx = fabs(x) - _scaleX*0.5;
	double dy = fabs(y) - _scaleY*0.5;

	if (dx > dy)
		x = nsgn(x)*0.5*_scaleX;
	else
		y = nsgn(y)*0.5*_scaleY;

	rotate(x, y, _theta);
	x += _posX;
	y += _posY;
}

void SolidBox::distanceNormal(double &nx, double &ny, double x, double y) const {
    x -= _posX;
    y -= _posY;
    rotate(x, y, -_theta);
    if (fabs(x) - _scaleX*0.5 > fabs(y) - _scaleY*0.5) {
        nx = nsgn(x);
        ny = 0.0;
    } else {
        nx = 0.0;
        ny = nsgn(y);
    }
    rotate(nx, ny, _theta);
}

double SolidSphere::distance(double x, double y) const {
    return length(x - _posX, y - _posY) - _scaleX*0.5;
}

void SolidSphere::closestSurfacePoint(double &x, double &y) const {
    globalToLocal(x, y);
    
	double r = length(x, y);
	if (r < 1e-4) {
		x = 0.5;
		y = 0.0;
	} else {
		x /= 2.0*r;
		y /= 2.0*r;
	}
    
	localToGlobal(x, y);
}

void SolidSphere::distanceNormal(double &nx, double &ny, double x, double y) const {
    x -= _posX;
    y -= _posY;
    float r = length(x, y);
    if (r < 1e-4) {
        nx = 1.0;
        ny = 0.0;
    } else {
        nx = x/r;
        ny = y/r;
    }
}
