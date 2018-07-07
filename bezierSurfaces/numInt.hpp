/*************************************************************************
* Author: Christopher Dubbs
* Email: dubbsc@oregonstate.edu
* Date: April 12, 2018
* Program Name: OpenMP: Numeric Integration
* Program #: 1
*************************************************************************/

#include <iostream>
#include <iomanip>
#include <omp.h>
#include <string>
#include <fstream>
  

#define REPITITIONS 5000

// TOP SURFACE
#define XMIN 0.0
#define XMAX 3.0
#define YMIN 0.0
#define YMAX 3.0

#define TOPZ00 0.0
#define TOPZ10 1.0
#define TOPZ20 0.0
#define TOPZ30 0.0

#define TOPZ01 1.0
#define TOPZ11 6.0
#define TOPZ21 1.0
#define TOPZ31 0.0

#define TOPZ02 0.0
#define TOPZ12 1.0
#define TOPZ22 0.0
#define TOPZ32 4.0

#define TOPZ03 3.0
#define TOPZ13 2.0
#define TOPZ23 3.0
#define TOPZ33 3.0

// BOTTOM SURFACE
#define BOTZ00 0.0
#define BOTZ10 -3.0
#define BOTZ20 0.0
#define BOTZ30 0.0

#define BOTZ01 -2.0
#define BOTZ11 10.0
#define BOTZ21 -2.0
#define BOTZ31 0.0

#define BOTZ02 0.0
#define BOTZ12 -5.0
#define BOTZ22 0.0
#define BOTZ32 -6.0

#define BOTZ03 -3.0
#define BOTZ13 2.0
#define BOTZ23 -8.0
#define BOTZ33 -3.0


// Function Prototypes 
double Height(int iu, int iv);











