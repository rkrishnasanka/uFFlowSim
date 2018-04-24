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
    
    unsigned char *image = new unsigned char[sizeX*sizeY*4];

    vector<SolidBody *> bodies;
    double x, y, w, h, d, u, v, t;
    double posX, posY, scaleX, scaleY, theta, velX, velY, velTheta;
    string cmd;

    // Accepts file from command line with parameters
    if (argc > 1) { 
        io::CSVReader<9> in(argv[1]);

        while (in.read_row(cmd, posX, posY, scaleX, scaleY, theta, velX, velY, velTheta)) {
            // Specifies fluid parameters
            if (cmd == "fluid") {
                x = posX;
                y = posY;
                w = scaleX;
                h = scaleY;
                d = theta;
                u = velX;
                v = velY;
                t = velTheta;
            } 
            // Specifies matrix properties
            else if (cmd == "matrix") {
                sizeX = (int)posX;
                sizeY = (int)posY;
                density = scaleX;
                timestep = scaleY;
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
        x = 0.45;
        y = 0.20;
        w = 0.15;
        h = 0.03;
        d = 1.0;
        u = 0.0;
        v = 3.0;

        bodies.push_back(new SolidSphere(0.6, 0.7, 0.4, 0.1, M_PI*0.25, 0.0, 0.0, 0.0));
        bodies.push_back(new SolidSphere(0.1, 0.2, 0.1, 0.1, M_PI*0.1, 0.0, 0.0, 0.0));
    }
    
    vector<const SolidBody *> cBodies;
    for (unsigned i = 0; i < bodies.size(); i++)
        cBodies.push_back(bodies[i]);

    FluidSolver * solver = new FluidSolver(sizeX, sizeY, density, cBodies);

    double time = 0.0;
    int iterations = 0;

    chrono::duration<double> difference_in_time;
   
    
    while (time < t) {

    	// Start time
        chrono::time_point<chrono::steady_clock> begin_time = chrono::steady_clock::now();

        // Main iteration loop
        for (int i = 0; i < 4; i++) {
            solver->addInflow(x, y, w, h, d, u, v);
            solver->update(timestep);
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

    // Print timing result
    printf("Total time: %.10f seconds.\n", difference_in_time.count());

    // Terminate
    return 0;
}
