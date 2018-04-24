#include "FluidSolver.h"

void FluidSolver::buildRhs() {
    double scale = 1.0/_hx;
    const uint8_t *cell = _d->cell();
    const uint8_t *body = _d->body();
    
    for (int y = 0, idx = 0; y < _h; y++) {
        for (int x = 0; x < _w; x++, idx++) {
            if (cell[idx] == CELL_FLUID) {
                _r[idx] = -scale*
                    (_u->volume(x + 1, y)*_u->at(x + 1, y) - _u->volume(x, y)*_u->at(x, y) +
                     _v->volume(x, y + 1)*_v->at(x, y + 1) - _v->volume(x, y)*_v->at(x, y));
                
                double vol = _d->volume(x, y);
                
                if (_bodies.empty())
                    continue;
                
                if (x > 0)
                    _r[idx] -= (_u->volume(x, y) - vol)*_bodies[body[idx -  1]]->velocityX(x*_hx, (y + 0.5)*_hx);
                if (y > 0)
                    _r[idx] -= (_v->volume(x, y) - vol)*_bodies[body[idx - _w]]->velocityY((x + 0.5)*_hx, y*_hx);
                if (x < _w - 1)
                    _r[idx] += (_u->volume(x + 1, y) - vol)*_bodies[body[idx +  1]]->velocityX((x + 1.0)*_hx, (y + 0.5)*_hx);
                if (y < _h - 1)
                    _r[idx] += (_v->volume(x, y + 1) - vol)*_bodies[body[idx + _w]]->velocityY((x + 0.5)*_hx, (y + 1.0)*_hx);
            } else
                _r[idx] = 0.0;
        }
    }
}

/* Entries of the pressure matrix are modified accordingly */
void FluidSolver::buildPressureMatrix(double timestep) {
    double scale = timestep/(_density*_hx*_hx);
    const uint8_t *cell = _d->cell();
    
    memset(_aDiag,  0, _w*_h*sizeof(double));
    memset(_aPlusX, 0, _w*_h*sizeof(double));
    memset(_aPlusY, 0, _w*_h*sizeof(double));

    for (int y = 0, idx = 0; y < _h; y++) {
        for (int x = 0; x < _w; x++, idx++) {
            if (cell[idx] != CELL_FLUID)
                continue;
            
            if (x < _w - 1 && cell[idx + 1] == CELL_FLUID) {
                double factor = scale*_u->volume(x + 1, y);
                _aDiag [idx    ] +=  factor;
                _aDiag [idx + 1] +=  factor;
                _aPlusX[idx    ]  = -factor;
            }
            if (y < _h - 1 && cell[idx + _w] == CELL_FLUID) {
                double factor = scale*_v->volume(x, y + 1);
                _aDiag [idx     ] +=  factor;
                _aDiag [idx + _w] +=  factor;
                _aPlusY[idx     ]  = -factor;
            }
        }
    }
}

void FluidSolver::buildPreconditioner() {
    const double tau = 0.97;
    const double sigma = 0.25;
    const uint8_t *cell = _d->cell();

    for (int y = 0, idx = 0; y < _h; y++) {
        for (int x = 0; x < _w; x++, idx++) {
            if (cell[idx] != CELL_FLUID)
                continue;
            
            double e = _aDiag[idx];

            if (x > 0 && cell[idx - 1] == CELL_FLUID) {
                double px = _aPlusX[idx - 1]*_precon[idx - 1];
                double py = _aPlusY[idx - 1]*_precon[idx - 1];
                e = e - (px*px + tau*px*py);
            }
            if (y > 0 && cell[idx - _w] == CELL_FLUID) {
                double px = _aPlusX[idx - _w]*_precon[idx - _w];
                double py = _aPlusY[idx - _w]*_precon[idx - _w];
                e = e - (py*py + tau*px*py);
            }

            if (e < sigma*_aDiag[idx])
                e = _aDiag[idx];

            _precon[idx] = 1.0/sqrt(e);
        }
    }
}

void FluidSolver::applyPreconditioner(double *dst, double *a) {
    const uint8_t *cell = _d->cell();
    
    for (int y = 0, idx = 0; y < _h; y++) {
        for (int x = 0; x < _w; x++, idx++) {
            if (cell[idx] != CELL_FLUID)
                continue;
            
            double t = a[idx];

            if (x > 0 && cell[idx -  1] == CELL_FLUID)
                t -= _aPlusX[idx -  1]*_precon[idx -  1]*dst[idx -  1];
            if (y > 0 && cell[idx - _w] == CELL_FLUID)
                t -= _aPlusY[idx - _w]*_precon[idx - _w]*dst[idx - _w];

            dst[idx] = t*_precon[idx];
        }
    }

    for (int y = _h - 1, idx = _w*_h - 1; y >= 0; y--) {
        for (int x = _w - 1; x >= 0; x--, idx--) {
            if (cell[idx] != CELL_FLUID)
                continue;
            
            double t = dst[idx];

            if (x < _w - 1 && cell[idx +  1] == CELL_FLUID)
                t -= _aPlusX[idx]*_precon[idx]*dst[idx +  1];
            if (y < _h - 1 && cell[idx + _w] == CELL_FLUID)
                t -= _aPlusY[idx]*_precon[idx]*dst[idx + _w];

            dst[idx] = t*_precon[idx];
        }
    }
}

double FluidSolver::dotProduct(double *a, double *b) {
    double result = 0.0;
    for (int i = 0; i < _w*_h; i++)
        result += a[i]*b[i];
    return result;
}

void FluidSolver::matrixVectorProduct(double *dst, double *b) {
    for (int y = 0, idx = 0; y < _h; y++) {
        for (int x = 0; x < _w; x++, idx++) {
            double t = _aDiag[idx]*b[idx];
            
            if (x > 0)
                t += _aPlusX[idx -  1]*b[idx -  1];
            if (y > 0)
                t += _aPlusY[idx - _w]*b[idx - _w];
            if (x < _w - 1)
                t += _aPlusX[idx]*b[idx +  1];
            if (y < _h - 1)
                t += _aPlusY[idx]*b[idx + _w];

            dst[idx] = t;
        }
    }
}

void FluidSolver::scaledAdd(double *dst, double *a, double *b, double s) {
    for (int i = 0; i < _w*_h; i++)
        dst[i] = a[i] + b[i]*s;
}

double FluidSolver::infinityNorm(double *a) {
    double maxA = 0.0;
    for (int i = 0; i < _w*_h; i++)
        maxA = max(maxA, fabs(a[i]));
    return maxA;
}

void FluidSolver::project(int limit) {
    memset(_p, 0,  _w*_h*sizeof(double));
    applyPreconditioner(_z, _r);
    memcpy(_s, _z, _w*_h*sizeof(double));
    
    double maxError = infinityNorm(_r);
    if (maxError < 1e-5)
        return;
    
    double sigma = dotProduct(_z, _r);
    
    for (int iter = 0; iter < limit; iter++) {
        matrixVectorProduct(_z, _s);
        double alpha = sigma/dotProduct(_z, _s);
        scaledAdd(_p, _p, _s, alpha);
        scaledAdd(_r, _r, _z, -alpha);
        
        maxError = infinityNorm(_r);
        if (maxError < 1e-5) {
            printf("Exiting solver after %d iterations, maximum error is %f\n", iter, maxError);
            return;
        }
        
        applyPreconditioner(_z, _r);
        
        double sigmaNew = dotProduct(_z, _r);
        scaledAdd(_s, _z, _s, sigmaNew/sigma);
        sigma = sigmaNew;
    }
    
    printf("Exceeded budget of %d iterations, maximum error was %f\n", limit, maxError);
}

void FluidSolver::applyPressure(double timestep) {
    double scale = timestep/(_density*_hx);
    const uint8_t *cell = _d->cell();
    
    for (int y = 0, idx = 0; y < _h; y++) {
        for (int x = 0; x < _w; x++, idx++) {
            if (cell[idx] != CELL_FLUID)
                continue;
            
            _u->at(x,     y    ) -= scale*_p[idx];
            _u->at(x + 1, y    ) += scale*_p[idx];
            _v->at(x,     y    ) -= scale*_p[idx];
            _v->at(x,     y + 1) += scale*_p[idx];
        }
    }
}

void FluidSolver::setBoundaryCondition() {
    const uint8_t *cell = _d->cell();
    const uint8_t *body = _d->body();
    
    for (int y = 0, idx = 0; y < _h; y++) {
        for (int x = 0; x < _w; x++, idx++) {
            if (cell[idx] == CELL_SOLID) {
                const SolidBody &b = *_bodies[body[idx]];
                
                _u->at(x, y) = b.velocityX(x*_hx, (y + 0.5)*_hx);
                _v->at(x, y) = b.velocityY((x + 0.5)*_hx, y*_hx);
                _u->at(x + 1, y) = b.velocityX((x + 1.0)*_hx, (y + 0.5)*_hx);
                _v->at(x, y + 1) = b.velocityY((x + 0.5)*_hx, (y + 1.0)*_hx);
            }
        }
    }
    
    for (int y = 0; y < _h; y++)
        _u->at(0, y) = _u->at(_w, y) = 0.0;
    for (int x = 0; x < _w; x++)
        _v->at(x, 0) = _v->at(x, _h) = 0.0;
}

void FluidSolver::update(double timestep) {
    _d->fillSolidFields(_bodies);
    _u->fillSolidFields(_bodies);
    _v->fillSolidFields(_bodies);
    
    setBoundaryCondition();
    
    buildRhs();
    buildPressureMatrix(timestep);
    buildPreconditioner();
    project(2000);
    applyPressure(timestep);
    
    _d->extrapolate();
    _u->extrapolate();
    _v->extrapolate();
    
    setBoundaryCondition();
    
    _d->advect(timestep, *_u, *_v, _bodies);
    _u->advect(timestep, *_u, *_v, _bodies);
    _v->advect(timestep, *_u, *_v, _bodies);
    
    _d->flip();
    _u->flip();
    _v->flip();
}

void FluidSolver::addInflow(double x, double y, double w, double h, double d, double u, double v) {
    _d->addInflow(x, y, x + w, y + h, d);
    _u->addInflow(x, y, x + w, y + h, u);
    _v->addInflow(x, y, x + w, y + h, v);
}

void FluidSolver::addOutflow(double x, double y, double w, double h, double d, double u, double v) {
    _d->addOutflow(x, y, x + w, y + h, d);
    _u->addOutflow(x, y, x + w, y + h, u);
    _v->addOutflow(x, y, x + w, y + h, v);
}

void FluidSolver::toImage(unsigned char *rgba) {
    for (int i = 0; i < _w*_h; i++) {
        /* Use fluid volume for nice anti aliasing */
        int shade = (int)((1.0 - _d->src()[i])*_d->volume(i % _w, i/_w)*255.0);
        shade = max(min(shade, 255), 0);
        
        rgba[i*4 + 0] = shade;
        rgba[i*4 + 1] = shade;
        rgba[i*4 + 2] = shade;
        rgba[i*4 + 3] = 0xFF;
    }
}
