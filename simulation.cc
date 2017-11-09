#include "simulation.h"
#include "constants.h"

// Constructor reads parameters from a config file.
Simulation::Simulation (void) {

  string tmp;
  string item;
  int pos;
  ifstream f_conf("in/parameters.txt");

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


// Simulate the decoding algorithm. Run simulation and save result.
void Simulation::run (void) {

  // Random number initialization
  rngInitialize();

  // Create output directory
  stringstream ssdir, ssdel;
  ssdir << fixed << setprecision(2) << "out/p=" << pMin << "-" << pMax << "_pu=" << pLeakUpMin << "-" << pLeakUpMax;
  ssdel << "exec rm -Rf " << ssdir.str().c_str();
  system (ssdel.str().c_str());

  system (ssdel.str().c_str());
  mkdir (ssdir.str().c_str(), S_IRWXU);

  // One repetition per lattice size
  for (int i = 0; i < (int) XSize.size(); i++) {
    int T = XSize[i]; 

    // The additional Z error is set to be the value between the min and max extremes -- the serial version of the simulator only uses a single value of Z
    double addZ = (double) (addlZMin + addlZMax) / 2.0;

    // One repetition per leakage parameter
    for (double pLeakUp = pLeakUpMax; pLeakUp >= pLeakUpMin - 1E-6; pLeakUp -= pLeakUpStep) {
      cout << "RUNNING SIMULATION:   XSize=" << XSize[i] << "   iterations=" << iterations << "  circuitErr=" << circuitErr << "   pMin=" << pMin << "   pMax=" << pMax << "  pStep=" << pStep << "  p2=" << p2Prob << "  q=" << qProb << "  pLeakUp=" << pLeakUp << "  pLeakDown=" << pLeakDown << "  addZ=" << addZ << "  corLeak=" << corLeak << endl;

      ofstream f_out;
      stringstream ssf;
      ssf << ssdir.str().c_str() << "/results_d=" << XSize[i] << "_pu=" << fixed << setprecision(2) << pLeakUp << ".txt";
      f_out.open (ssf.str().c_str());

      // One repetition per error probability
      for (double p = pMax; p >= pMin - 1E-6; p -= pStep) {

        // Initialize probabilities
        double p2 = p2Prob * p;           // probability of a two-qubit Pauli error
        double q = qProb * p;             // probability of measurement error
	double pu = pLeakUp * p;          // leakage probability
        double pd = pLeakDown * p;        // probabilty of return to computational state

        // Find edge weights by counting errors in a lattice with d = 3
        Counter ctr;
        ctr.initCount ();
	Lattice myLatticeCtr (3, 3, p, p2, q, pu, pd, 0, circuitErr, corLeak, &ctr);
	myLatticeCtr.simulate ();
        ctr.stopCount();

	// Iterate over fault paths to finalize the edge weights
        while (ctr.nextFaultPath()) {
          Lattice myLatticeErr(3, 3, p, p2, q, pu, pd, 0, circuitErr, corLeak, &ctr);
          myLatticeErr.simulate ();
	  bool success = myLatticeErr.success ();
	  if (!success) {
            ctr.printLastPath();
            cout << "Trace of uncorrectable error: " << endl << myLatticeErr.log.str() << endl;
	  }
          assert (success);
	  if (debug) cout << "SIMULATION REPORTS success=" << success << ". CURRENT LATTICE STATE:" << endl;
	  if (debug) myLatticeErr.printState();
	}

        // Initialize parameters for iterations
        int cntSucc = 0;

        // One repetition per iteration 
        for (int j = 1; j <= iterations; j++ ) {

          if (debug) cout << "NEW LATTICE WAS GENERATED:" << endl;
          ctr.initSilent();
          Lattice myLattice (XSize[i], T, p, p2, q, pu, pd, addZ, circuitErr, corLeak, &ctr); 
          if (debug) myLattice.printState();

          if (debug) cout << "SIMULATING CORRECTION WITH T=" << T << ":" << endl;
          myLattice.simulate ();

          bool success = myLattice.success();
          if (success) cntSucc++;
          if (debug) cout << "SIMULATION REPORTS success=" << success << ". CURRENT LATTICE STATE:" << endl;
          if (debug) myLattice.printState();
        }
        double psuc = (double)cntSucc/(double)iterations;
        double sig = std::sqrt(psuc*(1-psuc)/(double)(iterations-1));
        cout << fixed << setprecision(6) << fixed << "Success at error level p=" << 100 * p << "% is " << 100 * psuc << "% +- " << 100 * sig << "%." << endl;
        f_out << fixed << setprecision(6) << 100 * p << " " << setprecision(6) << 100 * psuc << " " << setprecision(6) << 100 * sig << endl;
      }
      f_out.close();
    }
  }
}
