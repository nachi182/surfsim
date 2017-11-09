/*
qecsim.cc
Martin Suchara
May 2014
*/

#include "simulation.h"
#include "constants.h"
#include <cstdlib>
#include <iostream>

using namespace std;


int main(int argc, char **argv) {

  if (argc != 1) {
    cout << "Parameters will be ignored." << endl;  
  }

  Simulation mySimulation;
  mySimulation.run ();

  return 0;
}
