#include "parseInput.h"

// Constructor reads parameters from a config file.
void ParseInput::read (char* infile) {

  string tmp;
  string item;
  int pos;
  ifstream f_conf(infile);

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  stringstream ss(tmp);
  while(ss >> item){
    XSize.push_back(atoi(item.c_str()));
  }

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  iterations = atoi(tmp.c_str());

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  circuitErr = (tmp.compare(0,3,"yes") == 0);
  if (tmp.compare(0,3,"yes") != 0 && tmp.compare(0,2,"no") != 0) {
    cout << "Error in configuration file. Expecting circuitErr = yes/no." << endl;
    exit(1);
  }

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  pMin = atof(tmp.c_str());

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  pMax = atof(tmp.c_str());

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  pStep = atof(tmp.c_str());

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  p2Prob = atof(tmp.c_str());

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  qProb = atof(tmp.c_str());

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  pLeakUpMin = atof(tmp.c_str());

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  pLeakUpMax = atof(tmp.c_str());

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  pLeakUpStep = atof(tmp.c_str());

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  pLeakDown = atof(tmp.c_str());

  getline(f_conf,tmp);
  pos = tmp.find(" = ");
  tmp = tmp.substr (pos + 3);
  if (tmp.compare(0,2,"no") == 0)
    corLeak = no;
  else if (tmp.compare(0,5,"quick") == 0)
    corLeak = quick;
  else if (tmp.compare(0,4,"gate") == 0)
    corLeak = gate;
  else if (tmp.compare(0,7,"circuit") == 0)
    corLeak = circuit;
  else {
    cout << "Error in configuration file. Expecting correctLeak = no/quick/gate/circuit." << endl;
    exit(1);
  }

  // Optional parameters -- additional Z type error on idle qubits
  addlZMin = 0;
  addlZMax = 0;
  addlZStep = 0;
  if (getline(f_conf,tmp)) {
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    addlZMin = atof(tmp.c_str());
  }

  if (getline(f_conf,tmp)) {
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    addlZMax = atof(tmp.c_str());
  }

  if (getline(f_conf,tmp)) {
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    addlZStep = atof(tmp.c_str());
  }

  f_conf.close();

  // Check for combinations of parameters that are not allowed
  if (pLeakUpMax > 1E-10 && !circuitErr) {
    cout << "Leakage can be only simulated in the circuit model." << endl;
    exit(1);
  }

}


// Print all combinations of parameters.
void ParseInput::print (void) {

  //cout << "#XSize circuitErr p p2Prob qProb pLeakUp pLeakDown corLeak" << endl;
  // One repetition per lattice size
  for (int i = 0; i < (int) XSize.size(); i++) {
    // One repetition per leakage parameter
    for (double pLeakUp = pLeakUpMax; pLeakUp >= pLeakUpMin - 1E-6; pLeakUp -= pLeakUpStep) {
      //double myConst = (log(0.007)-log(0.0002))/30;
      //for (double p = 0; p <= 30; p = p+1) {
      //  double pLog = exp(log(0.0002)+p*myConst);
      // One repetition per error probability
      for (double p = pMax; p >= pMin - 1E-6; p -= pStep) {
        if (addlZMin == 0 && addlZMax == 0 && addlZStep == 0) {
	  cout << XSize[i] << " " << circuitErr << " " << p << " " << p2Prob << " " << qProb << " " << pLeakUp << " " << pLeakDown << " " << corLeak << endl;
        } else {
          for (double addlZ = addlZMax; addlZ >= addlZMin - 1E-6; addlZ -= addlZStep)
	    cout << XSize[i] << " " << circuitErr << " " << p << " " << p2Prob << " " << qProb << " " << pLeakUp << " " << pLeakDown << " " << corLeak << " " << addlZ << endl;
        }
      }
    }
  }
}


int main(int argc, char **argv) {

  ParseInput myParse;
  myParse.read (argv[1]);
  myParse.print ();

  return 0;
}
