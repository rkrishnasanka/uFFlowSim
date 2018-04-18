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
*/

#define _USE_MATH_DEFINES

#include <algorithm>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <stack>

#include "./lodepng/lodepng.h"

#include "Functions.h"
#include "SolidBody.h"
#include "FluidQuantity.h"
#include "FluidSolver.h"

using namespace std;

int main(int argc, char * argv[]) {
    /* Play with these constants, if you want */
    const int sizeX = 128;
    const int sizeY = 128;
    
    const double density = 0.1;
    const double timestep = 0.005;
    
    unsigned char *image = new unsigned char[sizeX*sizeY*4];


    double x, y, w, h, d, u, v;
    if (argc > 1) {
        x = atof(argv[1]);
        y = atof(argv[2]);
        w = atof(argv[3]);
        h = atof(argv[4]);
        d = atof(argv[5]);
        u = atof(argv[6]);
        v = atof(argv[7]);
    } else {
        x = 0.0;
        y = 0.0;
        w = 1.0;
        h = 1.0;
        d = 1.0;
        u = 0.0;
        v = 2.0;
    }

    
    vector<SolidBody *> bodies;
    bodies.push_back(new SolidSphere(0.6, 0.7, 0.4, 0.1, M_PI*0.25, 0.0, 0.0, 0.0));
    bodies.push_back(new SolidSphere(0.1, 0.2, 0.1, 0.1, M_PI*0.1, 0.0, 0.0, 0.0));
    
    vector<const SolidBody *> cBodies;
    for (unsigned i = 0; i < bodies.size(); i++)
        cBodies.push_back(bodies[i]);

    FluidSolver *solver = new FluidSolver(sizeX, sizeY, density, cBodies);

    double time = 0.0;
    int iterations = 0;
    
    while (time < 2.0) {
        for (int i = 0; i < 4; i++) {
            // solver->addInflow(0.45, 0.2, 0.15, 0.03, 1.0, 0.0, 3.0);
            solver->addInflow(x, y, w, h, d, u, v);
            solver->update(timestep);
            time += timestep;
            fflush(stdout);
        }
        
        solver->toImage(image);
        
        char path[256];
        sprintf(path, "Frame%05d.png", iterations++);
        lodepng_encode32_file(path, image, sizeX, sizeY);
        
        for (unsigned i = 0; i < bodies.size(); i++)
            bodies[i]->update(timestep);
    }

    return 0;
}
