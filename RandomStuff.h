#ifndef RandomStuff_cpp
#define RandomStuff_cpp

//#define USE_MPI


#ifdef USE_MPI
#include <mpi.h>
#endif

#include "RngStream.h"


#ifdef USE_MPI
void  rngInitializeMPI();
#endif

#ifndef USE_MPI
void rngInitialize();
#endif

double myRandom();


void rngDelete();

#endif

