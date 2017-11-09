#ifndef LATTICE_H
#define LATTICE_H

#include <cstdlib>
#include <assert.h>
#include <vector>
#include <string>
#include <sstream>
#include <cfloat>
#include <math.h>
#include "operator.h"
#include "../blossom5-v2.05.src/PerfectMatching.h"
#include "../blossom5-v2.05.src/GEOM/GeomPerfectMatching.h"
#include "constants.h"
#include "RandomStuff.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace std;
using namespace boost;

// Typedef for a graph from the Bosst library
typedef adjacency_list<listS, vecS, undirectedS, no_property, property<edge_weight_t, double> > mygraph;


class Vertex {
 public:
  int x;
  int y;
  int t;
  int id;
  int syndromeID;

  Vertex (int newX, int newY, int newT, int newID, int newSyndromeID);
  ~Vertex (void) { }
};


class MyEdge {
 public:
  int    id1;
  int    id2;
  double dist;

  MyEdge (int newID1, int newID2, double newDist);
  ~MyEdge (void) { }
};


class Syndrome {
 public:
  bool    value;
  int     syndromeID;
  Vertex* pVert;
  Vertex* nearestUp;
  Vertex* nearestDown;
  Vertex* nearestLeft;
  Vertex* nearestRight;
  Vertex* nearestAfter;

  Syndrome  (int newID);
  ~Syndrome (void) { }
};


class Cube {
 public:
  int xSize;
  int ySize;
  int TSize;
  Syndrome ****syndromes;
  bool ***lkgHistory;
  int numSynd;
  double wa, wb, wc, wd, we, wf;
  corTech corLeak;

  Cube (int newXSize, int newYSize, int newT);
  ~Cube (void);
  void accumulateRelativeWeights (bool circErr, Counter* ctr);
  void finalizeRelativeWeights (double p1, double p2, double pm, bool circErr, Counter* ctr);
  void initNearest (void);
  void getVertices (vector<Vertex*> & verticesXXXX, vector<Vertex*> & verticesZZZZ);
  void getEdges (vector<Vertex*> & vertices, vector<MyEdge*> & edges);
  void getEdgesLkgDetBoost (vector<Vertex*> & vertices, vector<MyEdge*> & edges, pauli errType);
  void updateEdge(int v1, int v2, double wtAdd, mygraph & g);
  void getEdgesLkgDetFloydWarshall (vector<Vertex*> & vertices, vector<MyEdge*> & edges, pauli errType);
  double getDistance (int x1Loc, int y1Loc, int t1Loc, int x2Loc, int y2Loc, int t2Loc);
  double min_double (vector<double> & vd);
  int min_3int (int x, int y, int z);  
};


class Cell {
 public:
  Operator *ErrCurrent;

  Cell  (void);
  ~Cell (void);
};


class Lattice {
 public:
  int xSize;
  int ySize;
  double p;
  double p2, q;
  double pu, pd, addZ, zero;
  int T;
  bool circErr;
  corTech corLeak;
  Cell ***cells;
  ostringstream log;
  Counter *ctr;

  Lattice (int newXSize, int newT, double pError, double p2Error, double qError, double pLeakUp, double pLeakDown, double addlZ, bool circuitErr, corTech correctLeak, Counter *pCtr);
  ~Lattice (void);
  void simulate (void);
  bool phenomSyndromeXXXX (int xLoc, int yLoc, double qProb, bool noisy=true);
  bool phenomSyndromeZZZZ (int xLoc, int yLoc, double qProb, bool noisy=true);
  void lkgRedCircuit (Cell* data, double pProb, double puProb, double pdProb, bool *lkgFlagD);
  void circSyndromeXXXX_1 (int xLoc, int yLoc, double pProb, double puProb, double pdProb);
  void circSyndromeXXXX_2 (int xLoc, int yLoc, double pProb, double puProb, double pdProb);
  void circSyndromeXXXX_3 (int xLoc, int yLoc, double pProb, double puProb, double pdProb);
  void circSyndromeXXXX_4 (int xLoc, int yLoc, double pProb, double puProb, double pdProb);
  void circSyndromeXXXX_5 (int xLoc, int yLoc, double pProb, double puProb, double pdProb);
  bool circSyndromeXXXX_6 (int xLoc, int yLoc, double pProb, double qProb, double puProb, double pdProb, bool *lkgFlagA, bool *lkgFlagD);
  void circSyndromeZZZZ_1 (int xLoc, int yLoc, double pProb, double puProb, double pdProb);
  void circSyndromeZZZZ_2 (int xLoc, int yLoc, double pProb, double puProb, double pdProb);
  void circSyndromeZZZZ_3 (int xLoc, int yLoc, double pProb, double puProb, double pdProb);
  void circSyndromeZZZZ_4 (int xLoc, int yLoc, double pProb, double puProb, double pdProb);
  void circSyndromeZZZZ_5 (int xLoc, int yLoc, double pProb, double puProb, double pdProb);
  bool circSyndromeZZZZ_6 (int xLoc, int yLoc, double pProb, double qProb, double puProb, double pdProb, bool *lkgFlagA, bool *lkgFlagD);
  void callMatching (pauli errType, vector<Vertex*> & vertices, vector<MyEdge*> & edges);
  void correctLine (pauli errType, int x1Loc, int y1Loc, int x2Loc, int y2Loc);
  void correct (pauli errType, int xLoc, int yLoc);
  Operator* getLogical (string whichOp);
  bool success (void);
  void printState (void);
};


#endif
