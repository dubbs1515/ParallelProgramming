/******************************************************************************
* Author: Christopher Dubbs
* Date: May 15, 2018
* CS475
* Program #: 5
* Program Name: Vectorized Array Multiplication and Reduction
******************************************************************************/
#include <iostream>
#include <omp.h>
#include <cstdlib>
#include <cstdio>
#include <iomanip> 
#include <fstream> 
#include "simd.p5.h"


// Globals
float arrA[ARRSIZE];
float arrB[ARRSIZE];
float arrC[ARRSIZE];
double avgTimeNonSIMD, avgTimeSIMD, timeElapsed, peakTimeNonSIMD = 1000000.0, peakTimeSIMD=1000000.0; 

// Function Prototypes
void fillArray(float *array, int arrSize);
void basicMUl(float *arr1, float *arr2, float *arr3, float size);
float basicMulSum(float *arr1, float *arr2, int size);
void printResults();

/******************************************************************************
MAIN
******************************************************************************/
int main()
{

    // Seed rand
    unsigned seed = time(0); 
    srand(seed);        

    // Fill arrays
    fillArray(arrA, ARRSIZE);
    fillArray(arrB, ARRSIZE);          


// Check to ensure that OpenMP is defined on current system	
#ifndef _OPENMP
    fprintf(stderr, "OpenMP is not supported in this environment.\n");
	return 1;
#endif


    
// ARRAY MULTIPLICATION TESTS
#if TESTTYPE == 0 

    // NON-SSE 
    for(int i = 0; i < NUMTRIES; i++)
    {
        double timeStartNonSIMD = omp_get_wtime();     // Log start time

        basicMUl(arrA, arrB, arrC, ARRSIZE);

        double timeFinishNonSIMD = omp_get_wtime();	// Log finish time

        timeElapsed = timeFinishNonSIMD - timeStartNonSIMD; 

        // Track peak time NON-SIMD (i.e. fastest)
        if(timeElapsed < peakTimeNonSIMD)
        {
            peakTimeNonSIMD = timeElapsed;
        }

        avgTimeNonSIMD  += timeElapsed;
    }

    avgTimeNonSIMD /= NUMTRIES;    // Finish average time calculation

    // SSE
    for(int i = 0; i < NUMTRIES; i++)
    {
        double timeStartSIMD = omp_get_wtime();     // Log start time

        SimdMul(arrA, arrB, arrC, ARRSIZE);

        double timeFinishSIMD = omp_get_wtime();	// Log finish time

        timeElapsed = timeFinishSIMD - timeStartSIMD; 

        // Track peak time SIMD (i.e. fastest)
        if(timeElapsed < peakTimeSIMD)
        {
            peakTimeSIMD = timeElapsed; 
        }

        avgTimeSIMD += timeElapsed; 
    }

    avgTimeSIMD /= NUMTRIES; 

#endif  //end array multiplication tests

// REDUCTION TESTS
#if TESTTYPE == 1
    float throwAwaySum; 

    // NON-SSE
    for(int i = 0; i < NUMTRIES; i++)
    {
        double timeStartNonSIMD = omp_get_wtime();     // Log start time

        throwAwaySum = basicMulSum(arrA, arrB, ARRSIZE);
    
        double timeFinishNonSIMD = omp_get_wtime();	// Log finish time

        timeElapsed = timeFinishNonSIMD - timeStartNonSIMD; 

        // Track peak time NON-SIMD (i.e. fastest)
        if(timeElapsed < peakTimeNonSIMD)
        {
            peakTimeNonSIMD = timeElapsed; 
        }

        avgTimeNonSIMD += timeElapsed; 
    }

    avgTimeNonSIMD /= NUMTRIES;   // Finish average time calculation

    //std::cout << "\nNON-SSE Reduction: " << throwAwaySum << std::endl;    // For testing

    // SSE 
    for(int i = 0; i < NUMTRIES; i++)
    {
        double timeStartSIMD = omp_get_wtime();     // Log start time
       
       throwAwaySum = SimdMulSum(arrA, arrB, ARRSIZE);

        double timeFinishSIMD = omp_get_wtime();	// Log finish time

        timeElapsed = timeFinishSIMD - timeStartSIMD; 

        // Track peak time SIMD (i.e. fastest)
        if(timeElapsed < peakTimeSIMD)
        {
            peakTimeSIMD = timeElapsed; 
        }

        avgTimeSIMD += timeElapsed;
    }

    avgTimeSIMD /= NUMTRIES;     // Finish average time calculation

    //std::cout << "\n SSE Reduction: " << throwAwaySum << std::endl;   //For testing

#endif  // end reduction tests

    // Print results
    printResults(); 
    
    return 0;
}


/******************************************************************************
* Function Name: printResults
******************************************************************************/
void printResults()
{
    // Open file to print results
	std::ofstream results;
	results.open("results.txt", std::ios::out | std::ios::app);

    std::cout << std::fixed << std::setprecision(10);
    results << std::fixed  << std::setprecision(10);

    // Print Experiment Type and Headings (if applicable)
    if(ARRSIZE == 1024)
    {
        if(TESTTYPE == 0)
        {
            std::cout << "\n-----------------------------------------------------\n";
            results << "\n-----------------------------------------------------\n";
            std::cout << "|          ARRAY MULTIPLICATION EXPERIMENT          |" << std::endl; 
            results << "|          ARRAY MULTIPLICATION EXPERIMENT          |" << std::endl;
            std::cout << "-----------------------------------------------------\n";
            results << "-----------------------------------------------------\n";
        }
        else
        {
            std::cout << "\n------------------------------------------\n";
            results << "\n------------------------------------------\n";
            std::cout << "|          REDUCTION EXPERIMENT          |" << std::endl; 
            results << "|          REDUCTION EXPERIMENT          |" << std::endl;
            std::cout << "------------------------------------------\n";
            results << "------------------------------------------\n";
        }

        std::cout << std::left << std::setw(24) << "ARRAY SIZE ";
        results << std::left << std::setw(24) << "ARRAY SIZE ";
        std::cout << std::left << std::setw(24) << "SIMD AVG (sec)";
        results << std::left << std::setw(24) << "SIMD AVG (sec)";
        std::cout << std::left << std::setw(24) << "SIMD PEAK (sec)";
        results << std::left << std::setw(24) << "SIMD PEAK (sec)";
        std::cout << std::left << std::setw(24) << "NON-SIMD AVG (sec) ";
        results << std::left << std::setw(24) << "NON-SIMD AVG (sec) "; 
        std::cout << std::left << std::setw(24) << "NON-SIMD PEAK (sec) ";
        results << std::left << std::setw(24) << "NON-SIMD PEAK (sec) ";
        std::cout << std::left << std::setw(24) << "SPEED-UP ";
        results << std::left << std::setw(24) << "SPEED-UP ";
        std::cout << std::endl;
        results << std::endl;
    }

    // Print results
    std::cout << std::left << std::setw(24) << ARRSIZE;
    results << std::left << std::setw(24) << ARRSIZE; 
    std::cout << std::left << std::setw(24) << avgTimeSIMD;
    results << std::left << std::setw(24) << avgTimeSIMD;
    std::cout << std::left << std::setw(24) << peakTimeSIMD;
    results << std::left << std::setw(24) << peakTimeSIMD;
    std::cout << std::left << std::setw(24) << avgTimeNonSIMD;
    results << std::left << std::setw(24) << avgTimeNonSIMD;
    std::cout << std::left << std::setw(24) << peakTimeNonSIMD;
    results << std::left << std::setw(24) << peakTimeNonSIMD;
    std::cout << std::left << std::setw(24) << (peakTimeNonSIMD/peakTimeSIMD);
    results << std::left << std::setw(24) << (peakTimeNonSIMD/peakTimeSIMD);
    std::cout << std::endl;
    results << std::endl; 

    results.close();
}


/******************************************************************************
* Function Name:    basicMul
******************************************************************************/
void basicMUl(float *arr1, float *arr2, float *arr3, float size)
{
    for(int i=0; i < size; i++)
    {
        arr3[i] = arr1[i] * arr2[i];
    }
}


/******************************************************************************
* Function Name: basicMulSum
******************************************************************************/
float basicMulSum(float *arr1, float *arr2, int size)
{
    float sum = 0; 

    for(int i = 0; i < size; i++)
    {
        sum += arr1[i] * arr2[i];
    }

    return sum; 
}


/******************************************************************************
* Function Name: arrayFill
******************************************************************************/
void fillArray(float *array, int arrSize)
{
    for(int i = 0; i < arrSize; i++)
    {
        array[i] = (rand() % 100) + 1;      // random # b/n 1 and 100
    }
}