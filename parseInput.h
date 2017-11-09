#ifndef PARSEINPUT_H
#define PARSEINPUT_H

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
#include "colex.hh"
#include "constants.h"


using namespace std;


class ParseInput {
 public:
  vector<int> XSize;
  int iterations;
  bool circuitErr;
  double pMin, pMax, pStep;
  double p2Prob, qProb;
  double pLeakUpMin, pLeakUpMax, pLeakUpStep, pLeakDown;
  double addlZMin, addlZMax, addlZStep;
  corTech corLeak;

  ParseInput (void) { }
  ~ParseInput (void) { }

  void read (char* infile);
  void print (void);

};


#endif
