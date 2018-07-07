/*************************************************************************
* Author: Christopher Dubbs
* Email: dubbsc@oregonstate.edu
* Date: April 12, 2018
* Program Name: OpenMP: Numeric Integration
* Program #: 1
* Description: This program calculates the volume between two Bezier
	surfaces. 
*************************************************************************/
#include "numInt.hpp"

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

	double fullTileArea, megaComps, sumMegaComps, totVol;
	
	// Specify # of threads desired 
	omp_set_num_threads(NUMTHREADS);

	// Calculate the area for a full tile (middle tiles) (adjusted for edge and corners)
	fullTileArea = (((XMAX - XMIN) / (float)(NUMNODES-1)) * ((YMAX - YMIN) / (float)(NUMNODES-1)));
	
	// Open file to print results
	std::ofstream results;
	results.open("results.txt", std::ios::out | std::ios::app);

	// Specify # of threads desired 
	omp_set_num_threads(NUMTHREADS);

	double timeStart = omp_get_wtime();	// Log start

	// Parallelize using reduction to prevent partial sum addition problems bc of scheduler
	#pragma omp parallel for default(none), shared(fullTileArea), reduction(+:totVol)
	for(int i = 0; i < NUMNODES*NUMNODES; i++)
	{
		int ix, iy;
		double tileVol;

		// Determine cordinates
		ix = i % NUMNODES;
		iy = i / NUMNODES;

		// Caclulate tile volume
		tileVol = fullTileArea * Height(ix, iy);

		// Adjust volume for edge or corner if necessary
		if(ix == 0 || ix == NUMNODES-1)	// edge
		{
			tileVol /= 2.0; // Halve
		}

		if(iy == 0 || iy == NUMNODES-1)	// edge (accounts for corner as well)
		{
			tileVol /= 2.0; // Halve (or essentially quarter)
		}

		// Add to total volume
		totVol += tileVol;
	}
	// Back to single thread
	double timeFinish = omp_get_wtime();	// Log finish
	double executionTime = timeFinish - timeStart;
	megaComps = (double)(NUMNODES*NUMNODES) / (executionTime) / 1000000.0;

	// Output Results to Screen:
	std::cout << "Volume: " << std::fixed << std::setprecision(10) << totVol << std::endl; 
	std::cout << "MegaHeights Computed per Second: " << std::fixed << std::setprecision(10) << megaComps << std::endl;
	std::cout << "Execution Time: " << std::fixed << std::setprecision(10) << executionTime << std::endl;

	// Output Results to File:

	results << "\n# Threads: " << NUMTHREADS << std::endl; 
	results << "# Nodes (subdivisions): " << NUMNODES << std::endl;
	results << "Volume: " << std::fixed << std::setprecision(10) << totVol << std::endl; 
	results << "MegaHeights Computed per Second: " << std::fixed << std::setprecision(10) << megaComps << std::endl;
	results << "Execution Time: " << std::fixed << std::setprecision(10) << executionTime << std::endl;

	results.close();	// Close results file

	return 0;
}


/*************************************************************************
* Function Name: Height
* Description: Taken from program assignment notes.
*************************************************************************/
double Height(int iu, int iv)
{
	float u = (float)iu / (float)(NUMNODES-1);
	float v = (float)iv / (float)(NUMNODES-1);
	
	// Basis functions:
	float bu0 = (1.0-u) * (1.0-u) * (1.0-u);
	float bu1 = 3.0 * u * (1.0-u) * (1.0-u);
	float bu2 = 3.0 * u * u * (1.0-u);
	float bu3 = u * u * u;

	float bv0 = (1.0-v) * (1.0-v) * (1.0-v);
	float bv1 = 3.0 * v * (1.0-v) * (1.0-v);
	float bv2 = 3.0 * v * v * (1.0-v);
	float bv3 = v * v * v;

	float top =       bu0 * ( bv0*TOPZ00 + bv1*TOPZ01 + bv2*TOPZ02 + bv3*TOPZ03 )
                        + bu1 * ( bv0*TOPZ10 + bv1*TOPZ11 + bv2*TOPZ12 + bv3*TOPZ13 )
                        + bu2 * ( bv0*TOPZ20 + bv1*TOPZ21 + bv2*TOPZ22 + bv3*TOPZ23 )
                        + bu3 * ( bv0*TOPZ30 + bv1*TOPZ31 + bv2*TOPZ32 + bv3*TOPZ33 );

    float bot =       bu0 * ( bv0*BOTZ00 + bv1*BOTZ01 + bv2*BOTZ02 + bv3*BOTZ03 )
                        + bu1 * ( bv0*BOTZ10 + bv1*BOTZ11 + bv2*BOTZ12 + bv3*BOTZ13 )
                        + bu2 * ( bv0*BOTZ20 + bv1*BOTZ21 + bv2*BOTZ22 + bv3*BOTZ23 )
                        + bu3 * ( bv0*BOTZ30 + bv1*BOTZ31 + bv2*BOTZ32 + bv3*BOTZ33 );

    return top - bot;	// if the bottom surface sticks out above the top surface
						// then that contribution to the overall volume is negative			
}
