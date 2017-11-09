
#include "RandomStuff.h"


/// Global random stream
RngStream * rng;


#ifdef USE_MPI
void  rngInitializeMPI()
{ // Random number generator for a parallel machine
	const int rank = MPI::COMM_WORLD.Get_rank();
	int ss = time(NULL);
	MPI::COMM_WORLD.Bcast( & ss, 1, MPI::INT, 0 ); // to make sure seed is synchronized.
	unsigned long seed[6];
	seed[0] = ss % 12345678;
	seed[1] = seed[0];
	seed[2] = seed[0];
	seed[3] = seed[0];
	seed[4] = seed[0];
	seed[5] = seed[0];
	RngStream::SetPackageSeed( seed );
	rng = new RngStream;
	rng->IncreasedPrecis( true ); // full double-precision RNG
	for( int i = 0; i < 30; i++ ) // rank must be <= 2**30
		if( rank & (1<<i) )
			rng->AdvanceState( 127 + i, 0 );
}
#endif

#ifndef USE_MPI
void  rngInitialize()
{ // Random number generator for a serial machine
	unsigned long seed[6];
	seed[0] =  time(NULL) % 12345;
	seed[1] = seed[0];
	seed[2] = seed[0];
	seed[3] = seed[0];
	seed[4] = seed[0];
	seed[5] = seed[0];
	RngStream::SetPackageSeed( seed );
	rng = new RngStream;
}

#endif

double myRandom()
{
	return rng->RandU01();
}


void rngDelete()
{
 delete rng;
}

