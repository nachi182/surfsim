#ifndef SIMULATION_H
#define SIMULATION_H

#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "lattice.h"
#include "constants.h"


using namespace std;


class Simulation {
 public:
  vector<int> XSize;
  int iterations;
  bool circuitErr;
  double pMin, pMax, pStep;
  double p2Prob, qProb;
  double pLeakUpMin, pLeakUpMax, pLeakUpStep, pLeakDown;
  double addlZMin, addlZMax, addlZStep;
  corTech corLeak;

  Simulation (void);
  ~Simulation (void) { }
  void run (void);

};


#endif
