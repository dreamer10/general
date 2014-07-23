// Written by Peter Kutz.

#include "constants.h"

// TODO: Adjust the various settings below, and add more if you want.

const int theMillisecondsPerFrame = 10;


//const int theDim[3] = {2, 2, 1};

//#define _DEBUG

#ifdef _DEBUG
const int theDim[3] = {4, 4, 1};
#else
const int theDim[3] = {24, 24, 1};
//const int theDim[3] = {12, 12, 1};
#endif


const double theCellSize = 0.5;

