#ifndef stdafx_H	// make sure we dont define the same stuff repeatedly, otherwise get errors
#define stdafx_H

#pragma once

//#include <stdexcept>
#include <SDKDDKVer.h>
#include <stdio.h>
#include <cmath>
#include <iostream>  // I/O 
#include <fstream>   // file I/O
#include <iomanip>   // format manipulation
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <new>
#include <tchar.h>
//#include <amp.h>	// C++ AMP header file - only for VS2012, cant use this with Armadillo???
#include <ppl.h>

#define ARMA_64BIT_WORD
#include "armadillo"
using namespace arma;

using namespace Concurrency;   // goes with ppl.h

// define simple functions
#define isnan_agm(x) (x!=x)

// can define constants and functions here
double rot_angle = 0.61+0*3.142/2;	//DLD rotation angle to match trap axes
time_t seconds_start, seconds_end, seconds1, seconds2;

int max_group_time = 3400;	//maximum time between first and last event in group - should be 3400 bins of 25ps
//int dead_time = 400;		//time after 1 group to wait before looking for next group - should be 400 bins of 25ps
int T_sum = 3200;			// Is precisely the time taken for signal to travel 8cm, speed 1e6m/s, bins 25e-12 - should be 3200 bins of 25ps
//int tolerance = 200;		//tolerance in bins - should be 200 bins of 25ps

double bin_time = 25e-12;          // DLD bin size of 25 ps
int v_perp_x = 526320;				// signal propogation velocity m/s
int v_perp_y = 526320;

#endif