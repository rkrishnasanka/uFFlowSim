#include "FluidQuantity.h"

double FluidQuantity::lerp(double a, double b, double x) const {
    return a*(1.0 - x) + b*x;
}

double FluidQuantity::cerp(double a, double b, double c, double d, double x) const {
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

void FluidQuantity::rungeKutta3(double &x, double &y, double timestep, const FluidQuantity &u, const FluidQuantity &v) const {
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

void FluidQuantity::flip() {
    swap(_src, _dst);
}

const double * FluidQuantity::src() const {
    return _src;
}

const uint8_t * FluidQuantity::cell() const {
    return _cell;
}

const uint8_t * FluidQuantity::body() const {
    return _body;
}

double FluidQuantity::at(int x, int y) const {
    return _src[x + y*_w];
}

double FluidQuantity::volume(int x, int y) const {
    return _volume[x + y*_w];
}

double & FluidQuantity::at(int x, int y) {
    return _src[x + y*_w];
}

double FluidQuantity::lerp(double x, double y) const {
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

double FluidQuantity::cerp(double x, double y) const {
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

void FluidQuantity::backProject(double &x, double &y, const vector<const SolidBody *> &bodies) {
    int rx = min(max((int)(x - _ox), 0), _w - 1);
    int ry = min(max((int)(y - _oy), 0), _h - 1);
    
    if (_cell[rx + ry*_w] != CELL_FLUID) {
        x = (x - _ox)*_hx;
        y = (y - _oy)*_hx;
        bodies[_body[rx + ry*_w]]->closestSurfacePoint(x, y);
        x = x/_hx + _ox;
        y = y/_hx + _oy;
    }
}

void FluidQuantity::advect(double timestep, const FluidQuantity &u, const FluidQuantity &v,
        const vector<const SolidBody *> &bodies) {

    // Start time
    chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();
    
    #pragma omp parallel for
    for (int iy = 0, idx = 0; iy < _h; iy++) {
        for (int ix = 0; ix < _w; ix++, idx++) {
            if (_cell[idx] == CELL_FLUID) {
                double x = ix + _ox;
                double y = iy + _oy;
                
                rungeKutta3(x, y, timestep, u, v);
                
                backProject(x, y, bodies);
                
                _dst[idx] = cerp(x, y);
            }
        }
    }

    // End time
    chrono::time_point<chrono::steady_clock> end_time = chrono::steady_clock::now();

    // Time difference
    chrono::duration<double> difference_in_time = end_time - begin_time;

    // Print candidate time
    if (PRINT_CANDIDATES) printf("Candidate advect time: %.10f seconds.\n", difference_in_time.count());

    // Add to aggregate candidate time
    candidate_advect_time += difference_in_time;
}

void FluidQuantity::addInflow(double x0, double y0, double x1, double y1, double v) {
    /*
    Figure out indices of the matrix where the inlet conditions will be placed
    */
    int ix0 = (int)(x0/_hx - _ox);
    int iy0 = (int)(y0/_hx - _oy);
    int ix1 = (int)(x1/_hx - _ox);
    int iy1 = (int)(y1/_hx - _oy);

    // Start time
    chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();
    
    //Span though the idices and change the source matrix values
    #pragma omp parallel for
    for (int y = max(iy0, 0); y < min(iy1, _h); y++) {
        for (int x = max(ix0, 0); x < min(ix1, _h); x++) {

            //We calculate the length function to use in the cubic pulse
            //Length = sqrt(x^2 + y^2)
            double l = length(
                (2.0*(x + 0.5)*_hx - (x0 + x1))/(x1 - x0),
                (2.0*(y + 0.5)*_hx - (y0 + y1))/(y1 - y0)
            );

            //We use the length calculated above on the pulse here
            double vi = cubicPulse(l)*v;
            if (fabs(_src[x + y*_w]) < fabs(vi))
                _src[x + y*_w] = vi;
        }
    }

    // End time
    chrono::time_point<chrono::steady_clock> end_time = chrono::steady_clock::now();

    // Time difference
    chrono::duration<double> difference_in_time = end_time - begin_time;

    // Print candidate time
    if (PRINT_CANDIDATES) printf("Candidate addInFlow time: %.10f seconds.\n", difference_in_time.count());

    // Add to aggregate candidate time
    candidate_addInFlow_time += difference_in_time;
}

void FluidQuantity::addOutflow(double x0, double y0, double x1, double y1, double v) {
    int ix0 = (int)(x0/_hx - _ox);
    int iy0 = (int)(y0/_hx - _oy);
    int ix1 = (int)(x1/_hx - _ox);
    int iy1 = (int)(y1/_hx - _oy);
    
    for (int y = max(iy0, 0); y < min(iy1, _h); y++) {
        for (int x = max(ix0, 0); x < min(ix1, _h); x++) {
            double l = length(
                (2.0*(x + 0.5)*_hx - (x0 + x1))/(x1 - x0),
                (2.0*(y + 0.5)*_hx - (y0 + y1))/(y1 - y0)
            );
            double vi = cubicPulse(l)*v;
            if (fabs(_src[x + y*_w]) > fabs(vi))
                _src[x + y*_w] = vi;
        }
    }
}

void FluidQuantity::fillSolidFields(const vector<const SolidBody *> &bodies) {
    if (bodies.empty())
        return;

    // Start time
    chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();
    
    /* Compute distance field first */
    #pragma omp parallel for
    for (int iy = 0, idx = 0; iy <= _h; iy++) {
        for (int ix = 0; ix <= _w; ix++, idx++) {
            double x = (ix + _ox - 0.5)*_hx;
            double y = (iy + _oy - 0.5)*_hx;
            
            _phi[idx] = bodies[0]->distance(x, y);
            for (unsigned i = 1; i < bodies.size(); i++)
                _phi[idx] = min(_phi[idx], bodies[i]->distance(x, y));
        }
    }
    
    #pragma omp parallel for
    for (int iy = 0, idx = 0; iy < _h; iy++) {
        for (int ix = 0; ix < _w; ix++, idx++) {
            double x = (ix + _ox)*_hx;
            double y = (iy + _oy)*_hx;
            
            _body[idx] = 0;
            double d = bodies[0]->distance(x, y);
            for (unsigned i = 1; i < bodies.size(); i++) {
                double id = bodies[i]->distance(x, y);
                if (id < d) {
                    _body[idx] = i;
                    d = id;
                }
            }
            
            /* Compute cell volume from the four adjacent distance samples */
            int idxp = ix + iy*(_w + 1);
            _volume[idx] = 1.0 - occupancy(
                _phi[idxp],          _phi[idxp + 1],
                _phi[idxp + _w + 1], _phi[idxp + _w + 2]
            );
            
            /* Clamp dangerously small cell volumes - could break numerical
             * solver otherwise
             */
            if (_volume[idx] < 0.01)
                _volume[idx] = 0.0;
            
            bodies[_body[idx]]->distanceNormal(_normalX[idx], _normalY[idx], x, y);
            
            /* Solid cells are now defined as cells with zero fluid volume */
            if (_volume[idx] == 0.0)
                _cell[idx] = CELL_SOLID;
            else
                _cell[idx] = CELL_FLUID;
        }
    }

    // End time
    chrono::time_point<chrono::steady_clock> end_time = chrono::steady_clock::now();

    // Time difference
    chrono::duration<double> difference_in_time = end_time - begin_time;

    // Print candidate time
    if (PRINT_CANDIDATES) printf("Candidate fillSolidFields time: %.10f seconds.\n", difference_in_time.count());

    // Add to aggregate candidate time
    candidate_fillSolidFields_time += difference_in_time;
}

void FluidQuantity::fillSolidMask() {

    // Start time
    chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();

    #pragma omp parallel for
    for (int y = 1; y < _h - 1; y++) {
        for (int x = 1; x < _w - 1; x++) {
            int idx = x + y*_w;
            
            if (_cell[idx] == CELL_FLUID)
                continue;
            
            double nx = _normalX[idx];
            double ny = _normalY[idx];

            _mask[idx] = 0;
            if (nx != 0.0 && _cell[idx + sgn(nx)]    != CELL_FLUID)
                _mask[idx] |= 1;
            if (ny != 0.0 && _cell[idx + sgn(ny)*_w] != CELL_FLUID)
                _mask[idx] |= 2;
        }
    }

    // End time
    chrono::time_point<chrono::steady_clock> end_time = chrono::steady_clock::now();

    // Time difference
    chrono::duration<double> difference_in_time = end_time - begin_time;

    // Print candidate time
    if (PRINT_CANDIDATES) printf("Candidate fillSolidMask time: %.10f seconds.\n", difference_in_time.count());

    // Add to aggregate candidate time
    candidate_fillSolidMask_time += difference_in_time;
}

double FluidQuantity::extrapolateNormal(int idx) {
    double nx = _normalX[idx];
    double ny = _normalY[idx];
    
    double srcX = _src[idx + sgn(nx)];
    double srcY = _src[idx + sgn(ny)*_w];
    
    return (fabs(nx)*srcX + fabs(ny)*srcY)/(fabs(nx) + fabs(ny));
}

void FluidQuantity::freeNeighbour(int idx, stack<int> &border, int mask) {
    _mask[idx] &= ~mask;
    if (_cell[idx] != CELL_FLUID && _mask[idx] == 0)
        border.push(idx);
}

void FluidQuantity::extrapolate() {
    fillSolidMask();

    stack<int> border;

    // Start time
    chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();

    #pragma omp parallel for
    for (int y = 1; y < _h - 1; y++) {
        for (int x = 1; x < _w - 1; x++) {
            int idx = x + y*_w;

            if (_cell[idx] != CELL_FLUID && _mask[idx] == 0)
                border.push(idx);
        }
    }

    while (!border.empty()) {
        int idx = border.top();
        border.pop();

        _src[idx] = extrapolateNormal(idx);

        if (_normalX[idx - 1] > 0.0)
            freeNeighbour(idx -  1, border, 1);
        if (_normalX[idx + 1] < 0.0)
            freeNeighbour(idx +  1, border, 1);
        if (_normalY[idx - _w] > 0.0)
            freeNeighbour(idx - _w, border, 2);
        if (_normalY[idx + _w] < 0.0)
            freeNeighbour(idx + _w, border, 2);
    }

    // End time
    chrono::time_point<chrono::steady_clock> end_time = chrono::steady_clock::now();

    // Time difference
    chrono::duration<double> difference_in_time = end_time - begin_time;

    // Print candidate time
    if (PRINT_CANDIDATES) printf("Candidate extrapolate time: %.10f seconds.\n", difference_in_time.count());

    // Add to aggregate candidate time
    candidate_extrapolate_time += difference_in_time;
}
