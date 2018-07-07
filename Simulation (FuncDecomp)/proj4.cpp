/*************************************************************************
* Author: Christopher Dubbs
* email: dubbsc@oregonstate.edu
* Date: May 8, 2018
* Program #: 4
* Program Name: Functional Decomposition
*************************************************************************/
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <iomanip>


// Globals
int nowYear;                 // 2017 - 2022
int nowMonth;                // 0 - 11
float nowPrecip;             // inches of rain per month
float nowTemp;               // temp this month
float nowHeight;             // grain height in inches
float nowHeightNoAttack;     // to plot grain height without added agent effect
int nowManBearPigStrikes;    // Manbearpig eats farmhands to the detriment of next month's crop yield (grain growth)
int nowNumDeer;              // # of deer in current population
unsigned int seed;           // For ranf


const float GRAIN_GROWS_PER_MONTH = 8.0;    // (inches) 
const float ONE_DEER_EATS_PER_MONTH = 0.5;
const float MAN_BEAR_PIG_EFFECT = 0.2;      // Lost crop yield (graingrowth) due to a eaten farmhand

const float AVG_PRECIP_PER_MONTH = 6.0;     // average  (inches)
const float AMP_PRECIP_PER_MONTH = 6.0;     // plus or minus
const float RANDOM_PRECIP = 2.0;            // plus or minus noise

const float AVG_TEMP = 50.0;                // average  (degrees F)
const float AMP_TEMP = 20.0;                // plus or minus
const float RANDOM_TEMP = 10.0;             // plus or minus noise

const float MIDTEMP = 40.0;
const float MIDPRECIP = 10.0;   
            

// Function Prototypes
float Ranf(unsigned int *seedp, float low, float high);
int Ranf(unsigned int *seedp, int ilow, int ihigh);
void updateWeather();
void grainDeer();
void grain();
void manBearPig();          // Added agent
void watcher(); 
float tempFactorCalc();     // Helper for grain()... affects growth
float precipFactorCalc();   // Helper for grain()... affects growth
void stateOutput();         // Save output to file and print to terminal

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


    // Set Starting Values 
    nowMonth = 0;               // Start Month 
    nowYear = 2017;             // Start year
    nowNumDeer = 1;             // Initial Deer Population
    nowHeight = 1.0;            // (i.e. grain height inches)
    nowHeightNoAttack = 1.0;    // Used for an extra chart highlighting impact of added agent
    nowManBearPigStrikes = 1;   // Attacks for fictional previous month,  (affects month 0)
    seed = time(NULL);          // Seed for rand_r
    updateWeather();            // Initial weather 

    omp_set_num_threads(4);     // Specify desired # of threads (same as # sections)

    // Use OpenMP sections to divide work among threads
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            grainDeer();
        }
        #pragma omp section
        {
            grain();
        }
        #pragma omp section
        {
            watcher();
        }
        #pragma omp section
        {
            manBearPig();
        }
    }
    // Implied barrier of OpenMP section
    
    return 0; 
}


/*************************************************************************
* Function Name: grain
* Description: updates the height of the grain based on influentional 
    factors (e.g. temperature, precipitation, manBearPig Attacks)
*************************************************************************/
void grain()
{
    while(nowYear < 2023)
    {
        float tempHeight = nowHeight;                   
        float tempHeightNoAttack = nowHeight; 
        float newPrecipFactor = precipFactorCalc();     // Update precip factor
        float newTempFactor = tempFactorCalc();         // Update temp factor

        // Effect of temp and precip
        tempHeight += newTempFactor * newPrecipFactor * GRAIN_GROWS_PER_MONTH;

        // Effect of deer population
        tempHeight -= (float)nowNumDeer * ONE_DEER_EATS_PER_MONTH;
        tempHeightNoAttack = tempHeight;

        // Effect of Manbearpig attacks
        tempHeight -= ((float) nowManBearPigStrikes * MAN_BEAR_PIG_EFFECT); 
        
        // Ensure grain height is a non-negative value
        if(tempHeight < 0)
        {
            tempHeight = 0;
            tempHeightNoAttack = 0; // To keep the curves synchronized for the writeup 
        }

        // Barrier 1: DoneComputing barrier
        #pragma omp barrier 

        // Update current grain height
        nowHeight = tempHeight; 
        nowHeightNoAttack = tempHeightNoAttack;

        // Barrier 2: DoneAssigning barrier
        #pragma omp barrier

        //Barrier 3: DonePrinting barrier
        #pragma omp barrier
    }
}


/*************************************************************************
* Function Name:  grainDeer
* Description: This function updates the grainDeer population in 
    accordance with the height of the grain. 
*************************************************************************/
void grainDeer()
{
    while (nowYear < 2023)
    {
        // Calculate updated deer population
        float tempPop = nowNumDeer;   // Store for modification
        if((float) tempPop > nowHeight)
        {
            tempPop--;      // Increment deer population
        }
        else if((float) tempPop < nowHeight)
        {
            tempPop++;      // Decrement deer population
        }

        // Barrier 1: DoneComputing barrier
        #pragma omp barrier 

        // Update current deer population
        nowNumDeer = tempPop; 

        // Barrier 2: DoneAssigning barrier
        #pragma omp barrier

        //Barrier 3: DonePrinting barrier
        #pragma omp barrier
    }
}


/*************************************************************************
* Function Name: manBearPig
* Description: Manbearpig attacks local farmhands who tend the grain 
    fields. For each farmhand attacked per month, the crop yield (i.e. 
    grain growth) is decremented by 0.2 inches for the next month. 
*************************************************************************/
void manBearPig()
{
    while(nowYear < 2023)
    {   
        int high = 20, low = 1; 
        int tempMbpStrikes = Ranf(&seed, low, high);     // # of Manbearpig attacks

        // Barrier 1: DoneComputing barrier
        #pragma omp barrier 

        // Update Manbearpig attack data
        nowManBearPigStrikes = tempMbpStrikes;

        // Barrier 2: DoneAssigning barrier
        #pragma omp barrier

        //Barrier 3: DonePrinting barrier
        #pragma omp barrier
    }
}


/*************************************************************************
* Function Name: watcher
* Description: This function prints the current set of global state 
    variables, increments the month count, and uses the new month to 
    compute the new temperature and precipitation. 
*************************************************************************/
void watcher()
{
    while (nowYear < 2023)
    {
        // Barrier 1: DoneComputing barrier
        #pragma omp barrier 

        // Barrier 2: DoneAssigning barrier
        #pragma omp barrier

        // Output state
        stateOutput();

        // Increment month 
        nowMonth += 1;
        if (nowMonth == 12)   // Account for new year
        {
            nowMonth = 0; 
            nowYear++;
        }

        // Update weather (temperature & precipitation)
        updateWeather();

        //Barrier 3: DonePrinting barrier
        #pragma omp barrier
    }
}


/*************************************************************************
* Function Name: weatherUpdate
* Description: Used to update the temperature and precipitation factors.
* Citation: Taken from assignment notes. 
*************************************************************************/
void updateWeather()
{
    float ang = (30.0 * (float) nowMonth + 15.0) * (M_PI / 180.0);

    // Temperature
    float temp = AVG_TEMP - AMP_TEMP * cos(ang); 
    nowTemp = temp + Ranf(&seed, -RANDOM_TEMP, RANDOM_TEMP);

    // Precipitation
    float precip = AVG_PRECIP_PER_MONTH + AMP_PRECIP_PER_MONTH * sin(ang);
    nowPrecip = precip + Ranf(&seed, -RANDOM_PRECIP, RANDOM_PRECIP);
    if(nowPrecip < 0.0)
    {
        nowPrecip = 0.0;
    }
}


/*************************************************************************
* Function Name: Ranf (overloaded)
*************************************************************************/
float Ranf(unsigned int *seedp, float low, float high)
{
	float r = (float) rand_r(seedp);

    return (low + r * (high - low) / (float)RAND_MAX);
}

int Ranf(unsigned int *seedp, int ilow, int ihigh)
{
    float low = (float)ilow;
    float high = (float)ihigh + 0.9999f;
 
    return (int) (Ranf(seedp, low, high));
}


/*************************************************************************
* Function Name:  tempFactorCalc
* Description: Helper function to perform temperature factor calculation.
*************************************************************************/
float tempFactorCalc()
{
    return (std::exp(-(((nowTemp - MIDTEMP)/10.0) * ((nowTemp - MIDTEMP)/10.0))));
}


/*************************************************************************
* Function Name: precipFactorCalc
* Description: Helper function to perform the precipitation factor 
    calculation. 
*************************************************************************/
float precipFactorCalc()
{
    return (std::exp(-(((nowPrecip - MIDPRECIP)/10.0) * ((nowPrecip - MIDPRECIP)/10.0))));
}


/*************************************************************************
* Function Name: stateOutput
*************************************************************************/
void stateOutput()
{
    // Open file to print results
	std::ofstream states;
	states.open("states.txt", std::ios::out | std::ios::app);

    std::cout << std::fixed << std::setprecision(3);
    states << std::fixed  << std::setprecision(3);
    // Print Headings 
    if(nowYear == 2017 && nowMonth == 0)
    {
        std::cout << std::left << std::setw(10) << "YEAR ";
        states << std::left << std::setw(10) << "YEAR ";
        std::cout << std::left << std::setw(10) << "MONTH ";
        states << std::left << std::setw(10) << "MONTH ";
        std::cout << std::left << std::setw(24) << "PRECIPITATION (cm) ";
        states << std::left << std::setw(24) << "PRECIPITATION (cm) ";
        std::cout << std::left << std::setw(24) << "TEMPERATURE (C) ";
        states << std::left << std::setw(24) << "TEMPERATURE (C) "; 
        std::cout << std::left << std::setw(24) << "MANBEARPIGATTACKS ";
        states << std::left << std::setw(24) << "MANBEARPIGATTACKS ";
        std::cout << std::left << std::setw(24) << "HEIGHT W/O ATTACK";
        states << std::left << std::setw(24) << "HEIGHT W/O ATTACK ";
        std::cout << std::left << std::setw(24) << "GRAIN HEIGHT ";
        states << std::left << std::setw(24) << "GRAIN HEIGHT ";
        std::cout << std::left << std::setw(24) << "GRAIN DEER ";
        states << std::left << std::setw(24) << "GRAIN DEER ";
        std::cout << std::endl;
        states << std::endl;
    }

    std::cout << std::left << std::setw(10) << nowYear;
    states << std::left << std::setw(10) << nowYear;
    std::cout << std::left << std::setw(10) << nowMonth;
    states << std::left << std::setw(10) << nowMonth;
    std::cout << std::left << std::setw(24) << (nowPrecip * 2.54);
    states << std::left << std::setw(24) << (nowPrecip * 2.54);
    std::cout << std::left << std::setw(24) << ((5.0/9.0)*(nowTemp-32));
    states << std::left << std::setw(24) << ((5.0/9.0)*(nowTemp-32));
    std::cout << std::left << std::setw(24) << nowManBearPigStrikes;
    states << std::left << std::setw(24) << nowManBearPigStrikes;
    std::cout << std::left << std::setw(24) << nowHeightNoAttack;
    states << std::left << std::setw(24) << nowHeightNoAttack;
    std::cout << std::left << std::setw(24) << nowHeight;
    states << std::left << std::setw(24) << nowHeight;
    std::cout << std::left << std::setw(24) << nowNumDeer;
    states << std::left << std::setw(24) << nowNumDeer;
    std::cout << std::endl;
    states << std::endl;
  
    states.close();
}

