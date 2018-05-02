/*
Copyright (c) 2013 Benedikt Bitterli

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.
*/

#include <algorithm>
#include <math.h>
#include <cstring>
#include <stdio.h>
#include <chrono>
#include <omp.h>

using namespace std;

int loopIterations;

int NUM_THREADS = 1;

double length(double x, double y) {
    return sqrt(x*x + y*y);
}

double cubicPulse(double x) {
    x = min(fabs(x), 1.0);
    return 1.0 - x*x*(3.0 - 2.0*x);
}


class FluidQuantity {
    double *_src;
    double *_dst;

    int _w;
    int _h;
    double _ox;
    double _oy;
    double _hx;
    
    double lerp(double a, double b, double x) const {
        return a*(1.0 - x) + b*x;
    }

    double cerp(double a, double b, double c, double d, double x) const {
        double xsq = x*x;
        double xcu = xsq*x;
        
        double minV = min(a, min(b, min(c, d)));
        double maxV = max(a, max(b, max(c, d)));

        double t =
            a*(0.0 - 0.5*x + 1.0*xsq - 0.5*xcu) +
            b*(1.0 + 0.0*x - 2.5*xsq + 1.5*xcu) +
            c*(0.0 + 0.5*x + 2.0*xsq - 1.5*xcu) +
            d*(0.0 + 0.0*x - 0.5*xsq + 0.5*xcu);
        
        return min(max(t, minV), maxV);
    }
    
    void rungeKutta3(double &x, double &y, double timestep, const FluidQuantity &u, const FluidQuantity &v) const {
        double firstU = u.lerp(x, y)/_hx;
        double firstV = v.lerp(x, y)/_hx;

        double midX = x - 0.5*timestep*firstU;
        double midY = y - 0.5*timestep*firstV;

        double midU = u.lerp(midX, midY)/_hx;
        double midV = v.lerp(midX, midY)/_hx;

        double lastX = x - 0.75*timestep*midU;
        double lastY = y - 0.75*timestep*midV;

        double lastU = u.lerp(lastX, lastY);
        double lastV = v.lerp(lastX, lastY);
        
        x -= timestep*((2.0/9.0)*firstU + (3.0/9.0)*midU + (4.0/9.0)*lastU);
        y -= timestep*((2.0/9.0)*firstV + (3.0/9.0)*midV + (4.0/9.0)*lastV);
    }
    
public:
    FluidQuantity(int w, int h, double ox, double oy, double hx)
            : _w(w), _h(h), _ox(ox), _oy(oy), _hx(hx) {
        _src = new double[_w*_h];
        _dst = new double[_w*_h];
                
        memset(_src, 0, _w*_h*sizeof(double));
    }
    
    ~FluidQuantity() {
        delete[] _src;
        delete[] _dst;
    }
    
    void flip() {
        swap(_src, _dst);
    }
    
    const double *src() const {
        return _src;
    }
    
    double at(int x, int y) const {
        return _src[x + y*_w];
    }
    
    double &at(int x, int y) {
        return _src[x + y*_w];
    }
    
    double lerp(double x, double y) const {
        x = min(max(x - _ox, 0.0), _w - 1.001);
        y = min(max(y - _oy, 0.0), _h - 1.001);
        int ix = (int)x;
        int iy = (int)y;
        x -= ix;
        y -= iy;
        
        double x00 = at(ix + 0, iy + 0), x10 = at(ix + 1, iy + 0);
        double x01 = at(ix + 0, iy + 1), x11 = at(ix + 1, iy + 1);
        
        return lerp(lerp(x00, x10, x), lerp(x01, x11, x), y);
    }
    
    double cerp(double x, double y) const {
        x = min(max(x - _ox, 0.0), _w - 1.001);
        y = min(max(y - _oy, 0.0), _h - 1.001);
        int ix = (int)x;
        int iy = (int)y;
        x -= ix;
        y -= iy;
        
        int x0 = max(ix - 1, 0), x1 = ix, x2 = ix + 1, x3 = min(ix + 2, _w - 1);
        int y0 = max(iy - 1, 0), y1 = iy, y2 = iy + 1, y3 = min(iy + 2, _h - 1);
        
        double q0 = cerp(at(x0, y0), at(x1, y0), at(x2, y0), at(x3, y0), x);
        double q1 = cerp(at(x0, y1), at(x1, y1), at(x2, y1), at(x3, y1), x);
        double q2 = cerp(at(x0, y2), at(x1, y2), at(x2, y2), at(x3, y2), x);
        double q3 = cerp(at(x0, y3), at(x1, y3), at(x2, y3), at(x3, y3), x);
        
        return cerp(q0, q1, q2, q3, y);
    }
    
    void advect(double timestep, const FluidQuantity &u, const FluidQuantity &v) {
        int iy, ix, idx;


        for (iy = 0, idx = 0; iy < _h; iy++) {

            #pragma omp parallel for private(iy, idx, ix)
            for (ix = 0; ix < _w; ix++) {
                double x = ix + _ox;
                double y = iy + _oy;
                rungeKutta3(x, y, timestep, u, v);
                
                _dst[idx] = cerp(x, y);
                idx++;
            }
        }
    }
    
    void addInflow(double x0, double y0, double x1, double y1, double v) {
        int ix0 = (int)(x0/_hx - _ox);
        int iy0 = (int)(y0/_hx - _oy);
        int ix1 = (int)(x1/_hx - _ox);
        int iy1 = (int)(y1/_hx - _oy);
        int x, y;
        double l, vi;


        #pragma omp parallel for private(x, y, l, vi, ix0, iy0, ix1, iy1)
        for (y = max(iy0, 0); y < min(iy1, _h); y++) {
            for (x = max(ix0, 0); x < min(ix1, _h); x++) {
                l = length(
                    (2.0*(x + 0.5)*_hx - (x0 + x1))/(x1 - x0),
                    (2.0*(y + 0.5)*_hx - (y0 + y1))/(y1 - y0)
                );
                vi = cubicPulse(l)*v;
                if (fabs(_src[x + y*_w]) < fabs(vi))
                    _src[x + y*_w] = vi;
            }
        }
    }
};

class FluidSolver {
    FluidQuantity *_d;
    FluidQuantity *_u;
    FluidQuantity *_v;
    
    int _w;
    int _h;
    
    double _hx;
    double _density;
    
    double *_r;
    double *_p;
    
    void buildRhs() {
        int x, y, idx;
        double scale = 1.0/_hx;

        for (y = 0, idx = 0; y < _h; y++) {

            #pragma omp parallel for private(x, y, idx)
            for (x = 0; x < _w; x++) {
                _r[idx] = -scale*(_u->at(x + 1, y) - _u->at(x, y) +
                                  _v->at(x, y + 1) - _v->at(x, y));
                idx++;
            }
        }
    }
    
    void project(int limit, double timestep) {
        int x, y, iter, idx;
        double scale = timestep/(_density*_hx*_hx);
        
        double maxDelta;

        // omp_set_num_threads(NUM_THREADS);

    #pragma omp parallel
    {
        // if (omp_get_thread_num() == 0 && loopIterations == 0) 
        //     printf("Number of threads: %d\n", omp_get_num_threads());

        #pragma omp for private(y, idx, x)
        for (iter = 0; iter < limit; iter++) {
            maxDelta = 0.0;

            for (y = 0, idx = 0; y < _h; y++) {
                for (x = 0; x < _w; x++, idx++) {
                    idx = x + y*_w;
                    
                    double diag = 0.0, offDiag = 0.0;
                    
                    if (x > 0) {
                        diag    += scale;
                        offDiag -= scale*_p[idx - 1];
                    }
                    if (y > 0) {
                        diag    += scale;
                        offDiag -= scale*_p[idx - _w];
                    }
                    if (x < _w - 1) {
                        diag    += scale;
                        offDiag -= scale*_p[idx + 1];
                    }
                    if (y < _h - 1) {
                        diag    += scale;
                        offDiag -= scale*_p[idx + _w];
                    }

                    double newP = (_r[idx] - offDiag)/diag;
                    
                    maxDelta = max(maxDelta, fabs(_p[idx] - newP));
                    
                    _p[idx] = newP;
                }
            }

            if (maxDelta < 1e-4) {
                // printf("Exiting solver after %d iterations, maximum change is %f\n", iter, maxDelta);
                loopIterations += iter;
                iter = limit;
            }
        }
    }
        // printf("Exceeded budget of %d iterations, maximum change was %f\n", limit, maxDelta);
    }
    
    void applyPressure(double timestep) {
        int x, y, idx;
        idx = 0;

        double scale = timestep/(_density*_hx);
        
        // #pragma omp parallel for private(x, y)
        for (y = 0; y < _h; y++) {

            #pragma omp parallel for private(x, y, idx)
            for (x = 0; x < _w; x++) {
                _u->at(x,     y    ) -= scale*_p[idx];
                _u->at(x + 1, y    ) += scale*_p[idx];
                _v->at(x,     y    ) -= scale*_p[idx];
                _v->at(x,     y + 1) += scale*_p[idx];
                idx++;
            }
        }
        
        #pragma omp parallel for private(x, y)
        for (y = 0; y < _h; y++)
            _u->at(0, y) = _u->at(_w, y) = 0.0;

        #pragma omp parallel for private(x, y)
        for (x = 0; x < _w; x++)
            _v->at(x, 0) = _v->at(x, _h) = 0.0;
    }
    
public:
    FluidSolver(int w, int h, double density) : _w(w), _h(h), _density(density) {
        _hx = 1.0/min(w, h);
        
        _d = new FluidQuantity(_w,     _h,     0.5, 0.5, _hx);
        _u = new FluidQuantity(_w + 1, _h,     0.0, 0.5, _hx);
        _v = new FluidQuantity(_w,     _h + 1, 0.5, 0.0, _hx);
        
        _r = new double[_w*_h];
        _p = new double[_w*_h];
        
        memset(_p, 0, _w*_h*sizeof(double));
    }
    
    ~FluidSolver() {
        delete _d;
        delete _u;
        delete _v;
        
        delete[] _r;
        delete[] _p;
    }
    
    void update(double timestep) {
        buildRhs();
        project(600, timestep);
        applyPressure(timestep);
        
        _d->advect(timestep, *_u, *_v);
        _u->advect(timestep, *_u, *_v);
        _v->advect(timestep, *_u, *_v);
        
        _d->flip();
        _u->flip();
        _v->flip();
    }
    
    void addInflow(double x, double y, double w, double h, double d, double u, double v) {
        _d->addInflow(x, y, x + w, y + h, d);
        _u->addInflow(x, y, x + w, y + h, u);
        _v->addInflow(x, y, x + w, y + h, v);
    }
};

int main(int argc, char * argv[]) {
    const int sizeX = 256;
    const int sizeY = 256;
    
    const double density = 0.1;
    const double timestep = 0.005;

    if (argc > 1) NUM_THREADS = atoi(argv[1]);
    else NUM_THREADS = 1;

    omp_set_num_threads(NUM_THREADS);

    FluidSolver *solver = new FluidSolver(sizeX, sizeY, density);

    double time = 0.0;
    loopIterations = 0;
    chrono::duration<double> difference_in_time;
    
    while (time < 5.0) {
        chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();
        for (int i = 0; i < 4; i++) {
            solver->addInflow(0.45, 0.2, 0.15, 0.03, 1.0, 0.0, 3.0);
            solver->update(timestep);
            time += timestep;
            fflush(stdout);
        }
        chrono::time_point<chrono::steady_clock> end_time = chrono::steady_clock::now();
        difference_in_time += end_time - begin_time;
    }

    // Print timing results
    printf("%.10f\n", difference_in_time.count());
    // printf("Total iterations: %d\n", loopIterations);
    // printf("Average iterations per timestep: %d\n", (int)(loopIterations / (5.0 / timestep)));

    return 0;
}
