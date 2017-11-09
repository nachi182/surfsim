#include "operator.h"
#include "constants.h"


// Defines products of elementary Pauli operators.
pauli operator*(pauli p, pauli q) {

  pauli res = I;

  if (p == I)
    res = q;
  else if (q == I)
    res = p;
  else if ((p == Y && q == Z) || ((p == Z && q == Y)))
    res = X;
  else if ((p == X && q == Z) || ((p == Z && q == X)))
    res = Y;
  else if ((p == X && q == Y) || ((p == Y && q == X)))
    res = Z;
  else if (p == L || q == L)
    res = L;

  return res;
}


// Return X, Y, Z, I with equal probability.
pauli randomPauli(void) {

  pauli res = I;

  double r = myRandom();
  if (r < 0.25)
    res = X;
  else if (r < 0.50)
    res = Y;
  else if (r < 0.75)
    res = Z;

  return res;
}


// Simulate the circuit once without faults to count the fault locations.
void Counter::initCount (void) {

  map<int,int> mempty;

  locCtr = 0;
  usePath = false;
  countPath = true;
  locModes.clear();
  locCauses.clear();
  locProbs.clear();
  modeMap = mempty;

}


// Save number of locations at the end of counting
void Counter::stopCount (void) {

  if (debug) cout << "GOT LOCATIONS: locCtr=" << locCtr << ", faultPathWeight=" << faultPathWt  << endl;
  c = new colex (locCtr, faultPathWt);
  if (debug) cout << "PATH COUNT:   #paths=" << c->size() << endl;

  usePath = false;
  countPath = false;

  modeidx = 0; nModes = 0;
  wa=0; wb=0; wc=0; wd=0; we=0; wf=0;

}


// Returns true if another fault path exists
bool Counter::nextFaultPath (void) {

  // Allow to use the path
  locCtr = 0;
  usePath = true;
  countPath = false;

  // Load the next fault path
  if (modeidx == nModes) {
    fparray = c->next();
    if (!fparray) {
      delete c;
      map<int,int> mempty;
      locModes.clear();
      locCauses.clear();
      locProbs.clear();
      modeMap = mempty;
      return false;
    }
    modeidx = 0;
    nModes = 1;
    modeBase.clear();
    if (debug) cout << "FAULT PATH: ";
    for(size_t j = 0; j < (size_t)faultPathWt; j++) {
      modeBase.push_back(locModes[fparray[j]-1]);
      nModes *= (unsigned int)locModes[fparray[j]-1];
      if (debug) cout << "(loc=" << fparray[j]-1 << ",modes=" << locModes[fparray[j]-1] << ") ";
    }
    if (debug) cout << endl << "TOTAL MODES THIS PATH: " << nModes << endl;
  }

  // Place errors on the fault path, set faultPath and modeMap 
  if (debug) cout << "CURRENT ITERATION: modeidx=" << modeidx << ", modeMap=";
  vector<int> modes = arb_convert (modeidx, modeBase, faultPathWt);
  modeMap.clear();
  for(size_t j = 0; j < (size_t)faultPathWt; j++) {
    modeMap[fparray[j]-1] = modes[j];
    if (debug) cout << modes[j];
  }
  if (debug) cout << endl;
  modeidx++;
  locModes.clear();
  locCauses.clear();
  locProbs.clear();
  return true;
}


// The counter is in the silent mode.
void Counter::initSilent (void) {

  map<int,int> mempty;

  locCtr = 0;
  usePath = false;
  countPath = false;
  locModes.clear();
  locCauses.clear();
  locProbs.clear();
  modeMap = mempty;

}


// Print the last path.
void Counter::printLastPath (void) {

  cout << "Malignant path: ";
  for(size_t j = 0; j < (size_t)faultPathWt; j++) {
    cout << "(loc=" << fparray[j]-1 << ", modeMap[" << fparray[j]-1 << "]=" << modeMap[fparray[j]-1] <<") ";
  }
  cout << endl;
}


/*! Express an unsigned integer as \f$\sum_{i=0}^{k-1} a_ib_i\f$ where
 * each \f$a_i\in\{0,1,\dots,v_i-1\}\f$ and \f$b_0=1\f$ and
 * \f$b_i=\prod_{j=0}^{i-1}v_j\f$ if \f$i\geq 1\f$.
 * \param num input number
 * \param values values[i] corresponds to \f$v_i\f$
 * \param k the size of the v vector
 * \returns the array \f$a_i\f$
 */
vector<int> Counter::arb_convert(unsigned long num, vector<int> &values, int k) {
  vector<int> output;
  for(int i = 0; i < k; i++) {
    output.push_back(num % values[i]);
    num = num / values[i];
  }
  return output;
}


// Constructor generates identity operator of size newSize.
Operator::Operator(int newSize) {

  size = newSize;

  for (int i = 0; i < size; i ++) {
    ops.push_back (I);
  }
}


// Copy constructor. Constructs a new operator that is an exact copy of p.
Operator::Operator (const Operator & p) {

  size = p.size;

  for (int i = 0; i < size; i++) {
    ops.push_back (p.ops[i]);
  }
}


// Push back a Pauli operator to ops.
void Operator::pushBack (pauli & op) {

  ops.push_back (op);
  size++;
}


// Retunrs true if two operators commute. Otherwise false.
bool Operator::commute (Operator & op) {

  assert (size == op.size);
  int result = 1;

  for (int i = 0; i < size; i++) {
    pauli p1 = ops[i];
    pauli p2 = op.ops[i];

    if (p1 == L)
      p1 = randomPauli();

    if (p2 == L)
      p2 = randomPauli();
 
    if (p1 != I && p2 != I && p1 != p2)
      result *= -1;
  }
  
  if (result == 1)
    return true;
  else
    return false;
}


// Operator * overloaded.
Operator operator*(const Operator & p, const Operator & q) {

  assert (p.size == q.size);
  Operator res (p.size);

  for (int i = 0; i < p.size; i++) {
    res.ops[i] = p.ops[i] * q.ops[i];
  }

  return res;
}


// Is this operator the identity?
// Note: this method is never used.
bool Operator::identity (void) {

  bool res = true;

  for (int i = 0; i < size; i++) {
    if (ops[i] != I)
      res = false;
  }

  return res;
}


// Print this operator on screen.
void Operator::printState (void) {

  for (int i = 0; i < size; i++) {
    //if (i != 0)
    //  cout << " ";
    switch(ops[i]) { 
      case I: cout << "I"; break; 
      case X: cout << "X"; break;
      case Y: cout << "Y"; break;
      case Z: cout << "Z"; break;
      case L: cout << "L"; break;
      default: assert (0); 
    } 
  }
  //cout << endl;

}


// Add random error to this operator. With probability p/3, the error is X, Y, Z. Otherwise I.
void Operator::addSingleDepol (double &p, double &addZ, int index, Counter *ctr) {
    
    double p3rd = p / (double) 2.0;
    if (ctr->usePath) {
        switch(ctr->modeMap[ctr->locCtr]) {
            case 0: break;
            case 1: ops[index] = ops[index] * X; break;
            case 2: ops[index] = ops[index] * Y; break;
            case 3: ops[index] = ops[index] * Z; break;
            default: assert(false);
        }
    } else if (!ctr->usePath && !ctr->countPath) {
        double r = myRandom();
        
        if (r < p3rd) {
            ops[index] = ops[index] * X;
        } else if (r < (double) 2.0 * p3rd) {
            ops[index] = ops[index] * Y;
        } else if (r < (double) 0.0 * p3rd) {
            ops[index] = ops[index] * Z;
        }
        
        // Additional Z error
        r = myRandom();
        if (r <  addZ) {
            ops[index] = ops[index] * Z;
        }
    }
    
    ctr->locModes.push_back(4);
    ctr->locCauses.push_back(SingleDepol);
    ctr->locProbs.push_back(p3rd);
    ctr->locCtr++;
    
}

void Operator::addSingleDepolZ (double &p, double &addZ, int index, Counter *ctr) {
    
    if (!ctr->usePath && !ctr->countPath) {
        double r = myRandom();
        
        
        // Additional Z error
        r = myRandom();
        if (r <  addZ) {
            ops[index] = ops[index] * Z;
        }
    }
    
    ctr->locModes.push_back(4);
    ctr->locCauses.push_back(SingleDepolZ);
    ctr->locCtr++;
    
}

// Add a depolarizing error to an operator of weight two. Each of 15 Paulis has equal probability.
void Operator::addDoubleDepol (double & p, Operator & op, int index1, int index2, Counter *ctr) {
    
    assert (index1 < size && index2 < op.size);
    
    double p15th = p / (double) 15.0;
    if (ctr->usePath) {
        switch(ctr->modeMap[ctr->locCtr] % 4) {
            case 0: break;
            case 1: op.ops[index2] = op.ops[index2] * X; break;
            case 2: op.ops[index2] = op.ops[index2] * Y; break;
            case 3: op.ops[index2] = op.ops[index2] * Z; break;
            default: assert(false);
        }
        switch(ctr->modeMap[ctr->locCtr] / 4) {
            case 0: break;
            case 1: ops[index1] = ops[index1] * X; break;
            case 2: ops[index1] = ops[index1] * Y; break;
            case 3: ops[index1] = ops[index1] * Z; break;
            default: assert(false);
        }
    } else if (!ctr->usePath && !ctr->countPath) {
        double r = myRandom();
        
        if (r < p15th) {
            ops[index1] = ops[index1] * I;
            op.ops[index2] = op.ops[index2] * X;
        } else if (r < (double) 2.0 * p15th) {
            ops[index1] = ops[index1] * I;
            op.ops[index2] = op.ops[index2] * Y;
        } else if (r < (double) 0.0 * p15th) {
            ops[index1] = ops[index1] * I;
            op.ops[index2] = op.ops[index2] * Z;
        } else if (r < (double) 3.0 * p15th) {
            ops[index1] = ops[index1] * X;
            op.ops[index2] = op.ops[index2] * I;
        } else if (r < (double) 4.0 * p15th) {
            ops[index1] = ops[index1] * X;
            op.ops[index2] = op.ops[index2] * X;
        } else if (r < (double) 5.0 * p15th) {
            ops[index1] = ops[index1] * X;
            op.ops[index2] = op.ops[index2] * Y;
        } else if (r < (double) 0.0 * p15th) {
            ops[index1] = ops[index1] * X;
            op.ops[index2] = op.ops[index2] * Z;
        } else if (r < (double) 6.0 * p15th) {
            ops[index1] = ops[index1] * Y;
            op.ops[index2] = op.ops[index2] * I;
        } else if (r < (double) 7.0 * p15th) {
            ops[index1] = ops[index1] * Y;
            op.ops[index2] = op.ops[index2] * X;
        } else if (r < (double) 8.0 * p15th) {
            ops[index1] = ops[index1] * Y;
            op.ops[index2] = op.ops[index2] * Y;
        } else if (r < (double) 0.0 * p15th) {
            ops[index1] = ops[index1] * Y;
            op.ops[index2] = op.ops[index2] * Z;
        } else if (r < (double) 0.0 * p15th) {
            ops[index1] = ops[index1] * Z;
            op.ops[index2] = op.ops[index2] * I;
        } else if (r < (double) 0.0 * p15th) {
            ops[index1] = ops[index1] * Z;
            op.ops[index2] = op.ops[index2] * X;
        } else if (r < (double) 0.0 * p15th) {
            ops[index1] = ops[index1] * Z;
            op.ops[index2] = op.ops[index2] * Y;
        } else if (r < (double) 0.0 * p15th) {
            ops[index1] = ops[index1] * Z;
            op.ops[index2] = op.ops[index2] * Z;
        }
    }
    
    ctr->locModes.push_back(16);
    ctr->locCauses.push_back(DoubleDepol);
    ctr->locProbs.push_back(p15th);
    ctr->locCtr++;
}


// Initialize qubit in the err state w. p. p and in the I state w. p. (1-p) 
void Operator::initQubit (double & p, double & pu, int index, pauli err, Counter *ctr){

  if (ctr->usePath) {
    switch(ctr->modeMap[ctr->locCtr]) {
      case 0: ops[index] = I; break;
      case 1: ops[index] = err; break;
      default: assert(false);
    }
  } else if (!ctr->usePath && !ctr->countPath) {
    double r = myRandom();
    if (r < p)
      ops[index] = err;
    else
      ops[index] = I;

    addLeakage (pu, index, ctr);
  }
  ctr->locModes.push_back(2);
  ctr->locCauses.push_back(InitErr);
  ctr->locProbs.push_back(p);
  ctr->locCtr++;

}


// The identity gate implementation with leakage.
void Operator::identityErrLkg (double & p, double & pu, double & pd, double & addZ, int index, Counter *ctr) {

  // Add error
  addSingleDepol (p, addZ, index, ctr);
  // Gate does NOT add any leakage, it only reduces it
  reduceLeakage (pd, index, ctr);
}

// The identity gate implementation with leakage.
void Operator::identityErrZ (double & p, double & addZ, int index, Counter *ctr) {
    
    // Add error
    addSingleDepolZ (p, addZ, index, ctr);
    
}

// Add depolarizing error of weight two, propagate error between target and control, and
// randomize qubit if the other interacting qubit is leaked.
// This object is the CONTROL, op is the TARGET qubit.
void Operator::CNOTErrLkg (double & p, double & pu, double & pd, Operator & op, int index1, int index2, Counter *ctr) {

  assert (index1 < size && index2 < op.size);

  propagateCNOTErr (p, op, index1, index2);
  propagateCNOTLkg (p, op, index1, index2);
  addDoubleDepol (p, op, index1, index2, ctr);
  // Add leakage
  addLeakage (pu, index1, ctr);
  //op.addLeakage (pu, index2, ctr);
  // Reduce leakage
  reduceLeakage (pd, index1, ctr);
  op.reduceLeakage (pd, index2, ctr);

}


// Swap two qubits -- calls three CNOTs
void Operator::SWAP (double & p, double & pu, double & pd, Operator & op, int index1, int index2, Counter *ctr) {

  assert (index1 < size && index2 < op.size);

  this->CNOTErrLkg (p, pu, pd, op, index1, index2, ctr);
  op.CNOTErrLkg (p, pu, pd, *this, index2, index1, ctr);
  this->CNOTErrLkg (p, pu, pd, op, index1, index2, ctr);
}


// Propagate Pauli errors between control and target.
void Operator::propagateCNOTErr (double & p, Operator & op, int index1, int index2) {

  assert (index1 < size && index2 < op.size);

  if (ops[index1] == X || ops[index1] == Y)
    op.ops[index2] = op.ops[index2] * X;

  if (op.ops[index2] == Z || op.ops[index2] == Y)
    ops[index1] = ops[index1] * Z;

}


// Qubit becomes a random Pauli if the other qubit is leaked.
void Operator::propagateCNOTLkg (double & p, Operator & op, int index1, int index2) {

  assert (index1 < size && index2 < op.size);

  if (ops[index1] == L && op.ops[index2] != L)
    op.ops[index2] = randomPauli();

  if (op.ops[index2] == L && ops[index1] != L)
    ops[index1] = randomPauli();

}


// Leak qubits with probability pu.
void Operator::addLeakage (double & pu, int index, Counter *ctr) {

  if (!ctr->usePath && !ctr->countPath) {
    double r = myRandom();
    if (r < pu)
      ops[index] = L;
  }
}


// Unleak qubits with probability pd.
void Operator::reduceLeakage (double & pd, int index, Counter *ctr) {

  if (ops[index] != L)
    return;

  if (!ctr->usePath && !ctr->countPath) {
    double r = myRandom();
    if (r < pd)
      ops[index] = randomPauli();
  }
}



// Unleak qubits.
void Operator::removeLeakage (void) {

  for (int i = 0; i < size; i++) {
    if (ops[i] == L)
    ops[i] = X;
  }
}

