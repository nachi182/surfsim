#ifndef OPERATOR_H
#define OPERATOR_H

#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <vector>
#include <set>
#include <map>
#include "RandomStuff.h"
#include "colex.hh"

using namespace std;


enum pauli {I, X, Y, Z, L};
pauli operator*(pauli p, pauli q);
pauli randomPauli(void);

enum errCause {PhenSynd, EqLAnc, EqLData, MeasErr, SingleDepol, SingleDepolZ, DoubleDepol, InitErr, Lpu, PropLCNOT, MeasL};

class Counter {
 public:
  // Variables needed for fault location counting:
  static const int faultPathWt = 1;  // weight of the fault path
  int  locCtr;                // counter for fault locations
  bool usePath;               // use fault path or not?
  bool countPath;             // count error locations or not?
  vector<int> locModes;       // vector giving number of possible Pauli errors for corresponding fault location
  vector<errCause> locCauses; // vector giving the cause of the corresponding error
  vector<double> locProbs;    // vector giving the probability of the corresponding error
  map<int,int> modeMap;       // map from fault locations in faultPath to "mode", i.e. Pauli to apply

  unsigned long modeidx, nModes;
  colex *c;
  int *fparray;
  vector<int> modeBase;

  double wa, wb, wc, wd, we, wf; // edge weights for decoding are automatically generated by the counter

  Counter  (void) {};
  ~Counter (void) {};

  void initCount (void);
  void stopCount (void);
  bool nextFaultPath (void);
  void initSilent (void);
  void printLastPath (void);

 private:
  vector<int> arb_convert(unsigned long num, vector<int> &values, int k);
};

class Operator {
 public:
    int size;
    vector<pauli> ops;

    Operator (int newSize);
    Operator (const Operator & p);
    void pushBack (pauli & op);
    bool commute (Operator & op);
    friend Operator operator*(const Operator & p, const Operator & q);
    bool identity (void);
    void printState (void);
    void addSingleDepol (double & p, double & addZ, int index, Counter *ctr);
    void addSingleDepolZ (double & p, double & addZ, int index, Counter *ctr);
    void addDoubleDepol (double & p, Operator & op, int index1, int index2, Counter *ctr);
    void initQubit (double & p, double & pu, int index, pauli err, Counter *ctr);
    void identityErrLkg (double & p, double & pu, double & pd, double & addZ, int index, Counter *ctr);
    void identityErrZ (double & p, double & addZ, int index, Counter *ctr);
    void CNOTErrLkg (double & p, double & pu, double & pd, Operator & op, int index1, int index2, Counter *ctr);
    void SWAP (double & p, double & pu, double & pd, Operator & op, int index1, int index2, Counter *ctr);
    void propagateCNOTErr (double & p, Operator & op, int index1, int index2);
    void propagateCNOTLkg (double & p, Operator & op, int index1, int index2);
    void addLeakage (double & pu, int index, Counter *ctr);
    void reduceLeakage (double & pd, int index, Counter *ctr);
    void removeLeakage (void);
};


#endif