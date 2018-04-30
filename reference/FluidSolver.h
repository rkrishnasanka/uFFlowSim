#ifndef FLUIDSOLVER_H
#define FLUIDSOLVER_H

#define PRINT_CANDIDATES 0

#include <vector>
#include <chrono>
// #include <omp.h>

#include "Functions.h"
#include "SolidBody.h"
#include "FluidQuantity.h"

using namespace std;

class FluidSolver {
public:
    FluidQuantity *_d;
    FluidQuantity *_u;
    FluidQuantity *_v;
    
    int _w;
    int _h;
    
    double _hx;
    double _density;
    
    double *_r;
    double *_p;
    double *_z;
    double *_s;
    double *_precon;
    
    double *_aDiag;
    double *_aPlusX;
    double *_aPlusY;
    
    const vector<const SolidBody *> &_bodies;

    chrono::duration<double> candidate_buildRHS_time;
    chrono::duration<double> candidate_buildPressureMatrix_time;
    chrono::duration<double> candidate_buildPreconditioner_time;
    chrono::duration<double> candidate_applyPreconditioner_time;
    chrono::duration<double> candidate_dotProduct_time;
    chrono::duration<double> candidate_matrixVectorProduct_time;
    chrono::duration<double> candidate_scaleAdd_time;
    chrono::duration<double> candidate_infinityNorm_time;
    chrono::duration<double> candidate_applyPressure_time;
    chrono::duration<double> candidate_setBoundaryCondition_time;
    
    /* We now modify the right hand side to "blend" between solid and fluid
     * velocity based on the cell volume occupied by fluid.
     */
    void buildRhs();

    /* Entries of the pressure matrix are modified accordingly */
    void buildPressureMatrix(double timestep);
    void buildPreconditioner();
    void applyPreconditioner(double *dst, double *a);
    double dotProduct(double *a, double *b);
    void matrixVectorProduct(double *dst, double *b);
    void scaledAdd(double *dst, double *a, double *b, double s);
    double infinityNorm(double *a);
    void project(int limit);
    void applyPressure(double timestep);
    void setBoundaryCondition();
    
    /**
     * @brief Construct a new Fluid Solver object
     * 
     * @param w Width of the simulation area
     * @param h Height of the simulation area
     * @param density Density of the fluid
     * @param bodies List of solid bodies
     */
    FluidSolver(int w, int h, double density, const vector<const SolidBody *> &bodies)
                : _w(w), _h(h), _density(density), _bodies(bodies) {
        
        //From what it seems, this is the grid density 
        _hx = 1.0/min(w, h);
        
        _d = new FluidQuantity(_w,     _h,     0.5, 0.5, _hx);
        _u = new FluidQuantity(_w + 1, _h,     0.0, 0.5, _hx);
        _v = new FluidQuantity(_w,     _h + 1, 0.5, 0.0, _hx);
        
        _r = new double[_w*_h];
        _p = new double[_w*_h];
        _z = new double[_w*_h];
        _s = new double[_w*_h];
        _aDiag  = new double[_w*_h];
        _aPlusX = new double[_w*_h];
        _aPlusY = new double[_w*_h];
        _precon = new double[_w*_h];
    }
    
    ~FluidSolver() {
        delete _d;
        delete _u;
        delete _v;
        
        delete[] _r;
        delete[] _p;
        delete[] _z;
        delete[] _s;
        delete[] _aDiag;
        delete[] _aPlusX;
        delete[] _aPlusY;
        delete[] _precon;
    }
    
    void update(double timestep);
    void addInflow(double x, double y, double w, double h, double d, double u, double v);
    void addOutflow(double x, double y, double w, double h, double d, double u, double v);
    void toImage(unsigned char *rgba);
};

#endif
