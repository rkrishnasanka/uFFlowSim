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

Source for this code: https://github.com/tunabrain/incremental-fluids

This code also includes timing routines, from https://github.com/BU-EC-HPC-S16/
EC500-High-Performance-Computing/blob/master/ReferenceCode/n01TimingOpenMP/timi
ng_vector_mpi.cpp

*/

#define _USE_MATH_DEFINES

#include <algorithm>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <stack>
#include <chrono>

#include "./lodepng/lodepng.h"
#include "csv.h"

#include "Functions.h"
#include "SolidBody.h"
#include "FluidQuantity.h"
#include "FluidSolver.h"

using namespace std;

int main(int argc, char * argv[]) {
    /* Play with these constants, if you want */
    int sizeX = 128;
    int sizeY = 128;
    
    double density = 0.1;
    double timestep = 0.005;
    
    unsigned char * image = new unsigned char[sizeX*sizeY*4];

    vector<SolidBody *> bodies;
    double t;
    double xOut, yOut, wOut, hOut, dOut, uOut, vOut;
    double posX, posY, scaleX, scaleY, theta, velX, velY, velTheta;
    string cmd;

    // Accepts file from command line with parameters
    if (argc > 1) { 
        io::CSVReader<9> in(argv[1]);

        while (in.read_row(cmd, posX, posY, scaleX, scaleY, theta, velX, velY, velTheta)) {
            // Specifies source parameters
            if (cmd == "inflow") {
                xIn = posX;
                yIn = posY;
                wIn = scaleX;
                hIn = scaleY;
                dIn = theta;
                uIn = velX;
                vIn = velY;
            } 
            // Specifies sink parameters
            else if (cmd == "outflow") {
                xOut = posX;
                yOut = posY;
                wOut = scaleX;
                hOut = scaleY;
                dOut = theta;
                uOut = velX;
                vOut = velY;
            } 
            // Specifies matrix properties
            else if (cmd == "matrix") {
                sizeX = (int)posX;
                sizeY = (int)posY;
                density = scaleX;
                timestep = scaleY;
                t = theta;
            }
            // Adds sphere object to grid
            else if (cmd == "circle") {
                bodies.push_back(new SolidSphere(posX, posY, scaleX, scaleY, M_PI*theta, velX, velY, velTheta));
            } 
            // Adds box object to grid
            else if (cmd == "box") {
                bodies.push_back(new SolidBox(posX, posY, scaleX, scaleY, M_PI*theta, velX, velY, velTheta));                
            } 
            else {
                printf("Invalid command; ignoring line.\n");
            }
        }
    } else { // Default cases for debugging 
        bodies.push_back(new SolidSphere(0.6, 0.7, 0.4, 0.1, M_PI*0.25, 0.0, 0.0, 0.0));
        bodies.push_back(new SolidSphere(0.1, 0.2, 0.1, 0.1, M_PI*0.1, 0.0, 0.0, 0.0));
    }
    
    vector<const SolidBody *> cBodies;
    for (unsigned i = 0; i < bodies.size(); i++)
        cBodies.push_back(bodies[i]);

    FluidSolver * solver = new FluidSolver(sizeX, sizeY, density, cBodies);

    double time = 0.0;
    int iterations = 0;
    int loopIterations = 0;

    chrono::duration<double> difference_in_time;
   
    
    while (time < t) {

    	// Start time
        chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();


        // Main iteration loop
        for (int i = 0; i < 4; i++) {
            solver->addInflow(xIn, yIn, wIn, hIn, dIn, uIn, vIn);
            solver->update(timestep);

            loopIterations += solver->_iters;
            // printf("%d\n", loopIterations);

            time += timestep;
            fflush(stdout);
        }

        // End time
        chrono::time_point<chrono::steady_clock> end_time = chrono::steady_clock::now();
        
        // Add to total
        difference_in_time += end_time - begin_time;

        // Render image
        solver->toImage(image);
        char path[256];
        sprintf(path, "Frame%05d.png", iterations++);
        lodepng_encode32_file(path, image, sizeX, sizeY);
        
        for (unsigned i = 0; i < bodies.size(); i++)
            bodies[i]->update(timestep);
    }

    // Print timing results
    printf("Total time: %.10f seconds.\n", difference_in_time.count());
    printf("Total iterations: %d\n", loopIterations);
    printf("Average iterations per timestep: %d\n", (int)(loopIterations / (t / timestep)));

    // Print candidate results
    printf("Candidate time percentages (seconds/total):\n");

    printf("%-20f candidate_buildRHS_time\n", solver->candidate_buildRHS_time.count() / difference_in_time.count());
    printf("%-20f candidate_buildPressureMatrix_time\n", solver->candidate_buildPressureMatrix_time.count() / difference_in_time.count());
    printf("%-20f candidate_buildPreconditioner_time\n", solver->candidate_buildPreconditioner_time.count() / difference_in_time.count());
    printf("%-20f candidate_applyPreconditioner_time\n", solver->candidate_applyPreconditioner_time.count() / difference_in_time.count());
    printf("%-20f candidate_dotProduct_time\n", solver->candidate_dotProduct_time.count() / difference_in_time.count());
    printf("%-20f candidate_matrixVectorProduct_time\n", solver->candidate_matrixVectorProduct_time.count() / difference_in_time.count());
    printf("%-20f candidate_scaleAdd_time\n", solver->candidate_scaleAdd_time.count() / difference_in_time.count());
    printf("%-20f candidate_infinityNorm_time\n", solver->candidate_infinityNorm_time.count() / difference_in_time.count());
    printf("%-20f candidate_applyPressure_time\n", solver->candidate_applyPressure_time.count() / difference_in_time.count());
    printf("%-20f candidate_setBoundaryCondition_time\n", solver->candidate_setBoundaryCondition_time.count() / difference_in_time.count());

    printf("%-20f candidate_advect_time\n", (
        solver->_d->candidate_advect_time.count() + 
        solver->_u->candidate_advect_time.count() + 
        solver->_v->candidate_advect_time.count() 
        ) / difference_in_time.count()
    );

    printf("%-20f candidate_addInFlow_time\n", (
        solver->_d->candidate_addInFlow_time.count() + 
        solver->_u->candidate_addInFlow_time.count() + 
        solver->_v->candidate_addInFlow_time.count()
        ) / difference_in_time.count()
    );

    printf("%-20f candidate_fillSolidFields_time\n", (
        solver->_d->candidate_fillSolidFields_time.count() + 
        solver->_u->candidate_fillSolidFields_time.count() + 
        solver->_v->candidate_fillSolidFields_time.count()
        ) / difference_in_time.count()
    );

    printf("%-20f candidate_extrapolate_time\n", (
        solver->_d->candidate_extrapolate_time.count() + 
        solver->_u->candidate_extrapolate_time.count() + 
        solver->_v->candidate_extrapolate_time.count()
        ) / difference_in_time.count()
    );

    printf("%-20f candidate_fillSolidMask_time\n", (
        solver->_d->candidate_fillSolidMask_time.count() + 
        solver->_u->candidate_fillSolidMask_time.count() + 
        solver->_v->candidate_fillSolidMask_time.count()
        ) / difference_in_time.count()
    );

    printf("%-20f total\n", (
        solver->candidate_buildRHS_time.count() +
        solver->candidate_buildPressureMatrix_time.count() +
        solver->candidate_buildPreconditioner_time.count() +
        solver->candidate_applyPreconditioner_time.count() +
        solver->candidate_dotProduct_time.count() +
        solver->candidate_matrixVectorProduct_time.count() +
        solver->candidate_scaleAdd_time.count() +
        solver->candidate_infinityNorm_time.count() +
        solver->candidate_applyPressure_time.count() +
        solver->candidate_setBoundaryCondition_time.count() +
        solver->_d->candidate_advect_time.count() +
        solver->_d->candidate_addInFlow_time.count() +
        solver->_d->candidate_fillSolidFields_time.count() +
        solver->_d->candidate_extrapolate_time.count() +
        solver->_d->candidate_fillSolidMask_time.count() +
        solver->_u->candidate_advect_time.count() +
        solver->_u->candidate_addInFlow_time.count() +
        solver->_u->candidate_fillSolidFields_time.count() +
        solver->_u->candidate_extrapolate_time.count() +
        solver->_u->candidate_fillSolidMask_time.count() +
        solver->_v->candidate_advect_time.count() +
        solver->_v->candidate_addInFlow_time.count() +
        solver->_v->candidate_fillSolidFields_time.count() +
        solver->_v->candidate_extrapolate_time.count() +
        solver->_v->candidate_fillSolidMask_time.count()
        ) / difference_in_time.count()
    );

    // Terminate
    return 0;
}
