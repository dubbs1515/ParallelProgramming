/*************************************************************************
* Author: Christopher Dubbs
* Email: dubbsc@oregonstate.edu
* Date: April 12, 2018
* Program Name: OpenMP: False Sharing
* Program #: 3
*************************************************************************/
#include <iostream>
#include <iomanip>
#include <stdlib.h> 
#include <stdio.h> 
#include <omp.h>
#include <fstream>

struct s {
    float value;
    int pad[NUMPAD];    // Padding values will range from 0 to 16
}Array[4];


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

    unsigned int bigNum = 1000000000;

    omp_set_num_threads(NUMTHREADS);      // Specify # threads desired

    float temp;     // Used by fix 2

    double timeStart = omp_get_wtime();     // Log start time

    
#if FIXTYPE == 1 // Fix type 1: padding (Force values to new cache line, once NUMPAD is large enough)
    #pragma omp parallel for
    for (int i = 0; i < 4; i++)
    {
        for(int j = 0; j < bigNum; j++)
        {
            Array[i].value = Array[i].value + 2.0;      // Increment
        }
    }    
#endif  // Fix 1

#if FIXTYPE == 2 // Fix type 2: local variables (mem space on own stack)   
    #pragma omp parallel for default(none), shared (Array, bigNum), private(temp)
    for (int i = 0; i < 4; i++)
    {
        temp = Array[i].value;    // Use memory on own stack
        for(int j = 0; j < bigNum; j++)
        {
            temp = temp + 2.0;  // Increment  
        }
        Array[i].value = temp;  // Update array value 
    }
#endif // Fix 2

    double timeFinish = omp_get_wtime();    // Log finish time

	// Open file to print results
	std::ofstream results;
	results.open("results.txt", std::ios::out | std::ios::app);

    // Calculate performance
    double execTime = timeFinish - timeStart; 
    double megaIncrements = bigNum / execTime / 1000000.0;
    
    // Print Performance Metrics to Screen and Results File
    if(FIXTYPE == 1)
    {
        std::cout << std::endl << "Fix Type 1--Padding: " << std::endl;
        results << std::endl << "Fix Type 1--Padding: " << std::endl;
        std::cout << "Padding: " << NUMPAD << std::endl;
        results << "Padding: " << NUMPAD << std::endl; 
    }
    else
    {
        std::cout << std::endl << "Fix Type 2--Private Variable: " << std::endl;
        results << std::endl << "Fix Type 2--Private Variable: " << std::endl;
    }

    std::cout << "Number of Threads: " << NUMTHREADS << std::endl;
    results << "Number of Threads: " << NUMTHREADS << std::endl;
    std::cout << "Execution Time: " << std::fixed << std::setprecision(10) << execTime << std::endl;
    results << "Execution Time: " << std::fixed << std::setprecision(10) << execTime << std::endl;
    std::cout << "MegaIncrements Computed per Second: " << std::fixed << std::setprecision(10) << megaIncrements << std::endl << std::endl; 
    results << "MegaIncrements Computed per Second: " << std::fixed << std::setprecision(10) << megaIncrements << std::endl;

    results.close(); 

    return 0; 
}


