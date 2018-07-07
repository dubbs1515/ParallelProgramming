/*************************************************************************
* Author: Christopher Dubbs
* email: dubbsc@oregonstate.edu
* Date: April 22, 2018
* Program #: 2
* Program Name: N-body Problem (Coarse vs Fine & Static vs Dynamic)
*************************************************************************/
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <omp.h>

// Globals
const double G = 6.67300e-11;	// m^3 / (kg s^2)
const double EARTH_MASS = 5.9742e24; 	// kg
const double EARTH_DIAMETER = 12756000.32;	// meters
const double TIMESTEP = 1.0;	// seconds

#define NUMBODIES 100 
#define NUMSTEPS 200 

struct body {
	float mass;
	float x, y, z;	// position
	float vx, vy, vz;	// velocity
	float fx, fy, fz; 	// forces
	float xnew, ynew, znew;
	float vxnew, vynew, vznew;	
};

typedef struct body Body;

Body bodies[NUMBODIES];

// Function Prototypes 
float getDistanceSquared(Body *, Body *);
float getUnitVector(Body *, Body *, float *, float *, float *);
float ranf(float, float);
int ranf(int, int);


/*************************************************************************
MAIN
*************************************************************************/
int main()
{

// Check to ensure that OpenMP is defined on current system	
#ifndef _OPENMP
    fprintf(stderr, "OpenMP is not supported in this environment.\n");
	return 1;
#endif

    // Check for number of available processors
	omp_set_num_threads(NUMTHREADS);
	int numProc = omp_get_num_procs();
	std::cout << "Number of available processors: " << numProc << std::endl;

	
	for(int i =0; i < NUMBODIES; i++)
    {
        bodies[i].mass = EARTH_DIAMETER * ranf(0.5f, 10.0f);
        bodies[i].x = EARTH_DIAMETER * ranf(-100.0f, 100.0f);
        bodies[i].y = EARTH_DIAMETER * ranf(-100.0f, 100.0f);
        bodies[i].z = EARTH_DIAMETER * ranf(-100.0f, 100.0f);
        bodies[i].vx = ranf(-100.0f, 100.0f);
        bodies[i].vy = ranf(-100.0f, 100.0f);
        bodies[i].vz = ranf(-100.0f, 100.0f);
    }

    double timeStart = omp_get_wtime(); // Log start time

    for(int t = 0; t < NUMSTEPS; t++)
    {
		// Check to see if coarse and static scheduling are specified
		#if COARSE
		#if STATIC
			//std::cout << "Coarse + Static" << std::endl;
			#pragma omp parallel for default(none), shared(bodies), schedule(static)	// Coarse + Static
		#else 
			//std::cout << "Coarse + Dynamic" << std::endl;
			#pragma omp parallel for default(none), shared(bodies), schedule(dynamic)	// Coarse + Dynamic
		#endif	// Static
		#endif  // Coarse
    	for(int i = 0; i < NUMBODIES; i++)
    	{
    		float fx = 0.0;
			float fy = 0.0;
    		float fz = 0.0;
    		Body *bi = &bodies[i];

			// Check to see if Fine (i.e. !COARSE) and Static Scheduling are specified
			// Use reduction to prevent schedular complications for cumulative addition (if applicable)
			#if !COARSE
			#if STATIC
				//std::cout << "Fine + Static" << std::endl;
				#pragma omp parallel for default(none), shared(bodies, bi), reduction(+:fx,fy,fz), schedule(static)	// Fine + Static
			#else 
				//std::cout << "Fine + Dynamic" << std::endl;
				#pragma omp parallel for default(none), shared(bodies, bi), reduction(+:fx,fy,fz), schedule(dynamic)	// Fine + Dynamic
			#endif	// Static
			#endif  // Not Coarse (Fine)
    		for(int j = 0; j < NUMBODIES; j++)
    		{
    			if(j == 1) continue;

    			Body *bj = &bodies[j];

    			float rsqd = getDistanceSquared(bi, bj);

    			if(rsqd > 0.0)
    			{
					float f = G * bi->mass * bj->mass / rsqd; 
    				float ux, uy, uz;

    				getUnitVector(bi, bj, &ux, &uy, &uz);

    				fx += f * ux; 
    				fy += f * uy;
    				fz += f * uz;
    			}
    		}
    		float ax = fx / bodies[i].mass;
    		float ay = fy / bodies[i].mass;
    		float az = fz / bodies[i].mass;

    		bodies[i].xnew = bodies[i].x + bodies[i].vx * TIMESTEP + 0.5 * ax * TIMESTEP * TIMESTEP;
    		bodies[i].ynew = bodies[i].y + bodies[i].vy * TIMESTEP + 0.5 * ay * TIMESTEP * TIMESTEP;
    		bodies[i].znew = bodies[i].z + bodies[i].vz * TIMESTEP + 0.5 * az * TIMESTEP * TIMESTEP;

    		bodies[i].vxnew = bodies[i].vx + ax * TIMESTEP;
    		bodies[i].vynew = bodies[i].vy + ay * TIMESTEP;
    		bodies[i].vznew = bodies[i].vz + az * TIMESTEP; 
    	}
    	// Setup the state for the next animation step
    	for(int i = 0; i < NUMBODIES; i ++)
    	{
    		bodies[i].x = bodies[i].xnew;
    		bodies[i].y = bodies[i].ynew;
    		bodies[i].z = bodies[i].znew;
    		bodies[i].vx = bodies[i].vxnew;
    		bodies[i].vy = bodies[i].vynew;
    		bodies[i].vz = bodies[i].vznew;
    	}
    }		

    double timeFinish = omp_get_wtime();	// Log finish time
	 
    // Calculate Performance Metrics
	double execTime = timeFinish - timeStart;
	double megaBodies = (float) (NUMBODIES * NUMBODIES * NUMSTEPS) / (execTime)/1000000.0;


	// Open file to print results
	std::ofstream results;
	results.open("results.txt", std::ios::out | std::ios::app);


	// Print Performance Metrics to Screen and Results File
	#if COARSE	// Print Label
		std::cout << std::endl << "Coarse + ";
		results << std::endl << "Coarse + ";
	#else
		std::cout << std::endl << "Fine + ";
		results << std::endl << "Fine + ";
	#endif
	
	#if STATIC	// Print Label
		std::cout << "Static" << std::endl;
		results << "Static" << std::endl;
	#else
		std::cout << "Dynamic" << std::endl;
		results << "Dynamic" << std::endl;
	#endif

	std::cout << "Number of Threads: " << NUMTHREADS << std::endl; 
	results << "Number of Threads: " << NUMTHREADS << std::endl;

	std::cout << "Execution time: " << execTime << std::endl;
	results << "Execution time: " << execTime << std::endl;

	std::cout << "MegaBodies Computed per Second: " << megaBodies << std::endl << std::endl;
	results << "MegaBodies Computed per Second: " << megaBodies << std::endl << std::endl;
	
	results.close();

	return 0;
}


/*************************************************************************
* Function Name: getDistanceSquared
*************************************************************************/
float getDistanceSquared(Body *bi, Body *bj)
{
	float dx = bi->x - bj->x;
	float dy = bi->y - bj->y;
	float dz = bi->z - bj->z;

	return dx*dx + dy*dy + dz*dz;
}


/*************************************************************************
* Function Name: getUnitVector
*************************************************************************/
float getUnitVector(Body *from, Body *to, float *ux, float *uy, float *uz)
{
	float dx = to->x - from->x;
	float dy = to->y - from->y; 
	float dz = to->z - from->z; 

	float d = sqrt(dx*dx + dy*dy + dz*dz);

	if(d > 0.0)
	{
		dx /= d;
		dy /= d;
		dz /= d;
	}

	*ux = dx;
	*uy = dy;
	*uz = dz;

	return d; 
}


/*************************************************************************
* Function Name: ranf (float)
*************************************************************************/
float ranf(float low, float high)
{
	float r = (float) rand();	//0-RAND_MAX

	return (low + r * (high - low) / (float) RAND_MAX);
}


/*************************************************************************
* Function Name: ranf (int)
*************************************************************************/
int ranf(int ilow, int ihigh)
{
	float low = (float) ilow;
	float high = (float) ihigh + 0.9999f;

	return (int) (ranf(low, high));
}


