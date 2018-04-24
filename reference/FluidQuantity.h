#ifndef FLUIDQUANTITY_H
#define FLUIDQUANTITY_H

#define _USE_MATH_DEFINES

#include <stack>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <math.h>

#include "Functions.h"
#include "SolidBody.h"

using namespace std;

class FluidQuantity {
    double *_src;
    double *_dst;
    
    /* Distance field induced by solids.
     * Since this is used to compute the cell volumes, the samples are offset
     * by (-0.5, -0.5) from the samples in _src and the grid is one larger in
     * each dimension. This way, each sample of fluid quantity has four samples
     * of the distance function surrounding it - perfect for computing the
     * cell volumes.
     */
    double *_phi;
    /* Fractional cell volume occupied by fluid */
    double *_volume;
    
    double *_normalX;
    double *_normalY;
    uint8_t *_cell;
    uint8_t *_body;
    uint8_t *_mask;

    int _w;
    int _h;
    double _ox;
    double _oy;
    double _hx;
    
    double lerp(double a, double b, double x) const;
    double cerp(double a, double b, double c, double d, double x) const;
    void rungeKutta3(double &x, double &y, double timestep, const FluidQuantity &u, const FluidQuantity &v) const;
    
public:
    FluidQuantity(int w, int h, double ox, double oy, double hx)
            : _w(w), _h(h), _ox(ox), _oy(oy), _hx(hx) {
        _src = new double[_w*_h];
        _dst = new double[_w*_h];
        
        /* Make distance grid one larger in each dimension */
        _phi = new double[(_w + 1)*(_h + 1)];
        _volume  = new double[_w*_h];
        _normalX = new double[_w*_h];
        _normalY = new double[_w*_h];
        
        _cell = new uint8_t[_w*_h];
        _body = new uint8_t[_w*_h];
        _mask = new uint8_t[_w*_h];
        
        for (int i = 0; i < _w*_h; i++) {
            _cell[i] = CELL_FLUID;
            _volume[i] = 1.0;
        }
        
        memset(_src, 0, _w*_h*sizeof(double));
    }
    
    ~FluidQuantity() {
        delete[] _src;
        delete[] _dst;
        
        delete[] _phi;
        delete[] _volume;
        delete[] _normalX;
        delete[] _normalY;
        
        delete[] _cell;
        delete[] _body;
        delete[] _mask;
    }
    
    void flip();
    const double *src() const;
    const uint8_t *cell() const;
    const uint8_t *body() const;
    double at(int x, int y) const;
    double volume(int x, int y) const;
    double &at(int x, int y);
    double lerp(double x, double y) const;
    double cerp(double x, double y) const;
    void backProject(double &x, double &y, const vector<const SolidBody *> &bodies);
    void advect(double timestep, const FluidQuantity &u, const FluidQuantity &v,
            const vector<const SolidBody *> &bodies);
    void addInflow(double x0, double y0, double x1, double y1, double v);
    void addOutflow(double x0, double y0, double x1, double y1, double v);
    void fillSolidFields(const vector<const SolidBody *> &bodies);
    void fillSolidMask();
    double extrapolateNormal(int idx);
    void freeNeighbour(int idx, stack<int> &border, int mask);   
    void extrapolate();
};

#endif
