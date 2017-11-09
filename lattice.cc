#include "lattice.h"
#include "constants.h"


// Constructor of a vertex in a graph used in minimum weight matching.                
Vertex::Vertex (int newX, int newY, int newT, int newID, int newSyndromeID) {
  x  = newX;
  y  = newY;
  t  = newT;
  id = newID;
  syndromeID = newSyndromeID;
}


// Constructor of an edge: id1 and id2 are indexes of vertices and dist is their distance.
MyEdge::MyEdge (int newID1, int newID2, double newDist) {
  id1  = newID1;
  id2  = newID2;
  dist = newDist;
}


// Constructor generates a single syndrome residing in the 3-D cube.
Syndrome::Syndrome (int newID) {
  value = false;
  syndromeID = newID;
  pVert = NULL;
}


// Constructor generates the 3-D cube.
Cube::Cube(int newXSize, int newYSize, int newT) {

  xSize = newXSize;
  ySize = newYSize;
  TSize = newT;

  // Initialize syndrome history data structures
  int newXID = 0;
  int newZID = 0;
  syndromes = new Syndrome***[xSize];
  for (int i = 0; i < xSize; i++) {
    syndromes[i] = new Syndrome**[ySize];
    for (int j = 0; j < ySize; j++) {
      if ((i%2 && !(j%2)) || (!(i%2) && j%2))
        continue;
      syndromes[i][j] = new Syndrome*[TSize];
      for (int t = 0; t < TSize; t++) {
        if (i%2 && j%2)
          syndromes[i][j][t] = new Syndrome (newZID++);
        else
          syndromes[i][j][t] = new Syndrome (newXID++);
      }
    }
  }
  assert (newXID == newZID);
  numSynd = newXID;

  // Initialize leakage history
  lkgHistory = new bool**[xSize];
  for (int i = 0; i < xSize; i++) {
    lkgHistory[i] = new bool*[ySize];
    for (int j = 0; j < ySize; j++) {
      lkgHistory[i][j] = new bool[TSize];
      for (int t = 0; t < TSize; t++) {
        lkgHistory[i][j][t] = false;
      }
    }
  }
}


// Destructor of the 3-D Cube.
Cube::~Cube(void) {

  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      if ((i%2 && !(j%2)) || (!(i%2) && j%2))
        continue;
      for (int t = 0; t < TSize; t++) {
        delete syndromes[i][j][t];
      }
      delete[] syndromes[i][j];
    }
    delete[] syndromes[i];
  }
  delete[] syndromes; 

  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      delete[] lkgHistory[i][j];
    }
    delete[] lkgHistory[i];
  }
  delete[] lkgHistory; 

}


// Obtain weights of the six types of edges. Labeling of edges a, b, c, d, and e is consistent with Fig. 4 in arXiv:1009.3686v1.
void Cube::accumulateRelativeWeights (bool circErr, Counter* ctr) {

  assert (ctr->usePath && !ctr->countPath);
  assert (ctr->faultPathWt == 1);  

  // Update wa, wb, wc, etc.
  if (syndromes[4][4][1]->value && !syndromes[4][4][2]->value) {
    ctr->wa += ctr->locProbs[ctr->fparray[0]-1];
    if (debug) cout << "Adding to wa: p=" << ctr->locProbs[ctr->fparray[0]-1] << " cause=" << ctr->locCauses[ctr->fparray[0]-1] << " loc=" << ctr->fparray[0]-1 << endl;
  }
  if (syndromes[4][2][1]->value && !syndromes[4][2][0]->value && syndromes[4][4][1]->value && !syndromes[4][4][0]->value) {
    ctr->wb += ctr->locProbs[ctr->fparray[0]-1];
    if (debug) cout << "Adding to wb: p=" << ctr->locProbs[ctr->fparray[0]-1] << " cause=" << ctr->locCauses[ctr->fparray[0]-1] << " loc=" << ctr->fparray[0]-1 << endl;
  }
  if (syndromes[4][2][1]->value && !syndromes[4][2][0]->value && syndromes[4][4][2]->value && !syndromes[4][4][1]->value) {
    ctr->wc += ctr->locProbs[ctr->fparray[0]-1];
    if (debug) cout << "Adding to wc: p=" << ctr->locProbs[ctr->fparray[0]-1] << " cause=" << ctr->locCauses[ctr->fparray[0]-1] << " loc=" << ctr->fparray[0]-1 << endl;
  }
  if (syndromes[2][4][1]->value && !syndromes[2][4][0]->value && syndromes[4][4][1]->value && !syndromes[4][4][0]->value) {
    ctr->wd += ctr->locProbs[ctr->fparray[0]-1];
    if (debug) cout << "Adding to wd: p=" << ctr->locProbs[ctr->fparray[0]-1] << " cause=" << ctr->locCauses[ctr->fparray[0]-1] << " loc=" << ctr->fparray[0]-1 << endl;
  }
  if (syndromes[2][4][1]->value && !syndromes[2][4][0]->value && syndromes[4][4][2]->value && !syndromes[4][4][1]->value) {
    ctr->we += ctr->locProbs[ctr->fparray[0]-1];
    if (debug) cout << "Adding to we: p=" << ctr->locProbs[ctr->fparray[0]-1] << " cause=" << ctr->locCauses[ctr->fparray[0]-1] << " loc=" << ctr->fparray[0]-1 << endl;
  }
  if (syndromes[4][2][1]->value && !syndromes[4][2][0]->value && syndromes[2][4][2]->value && !syndromes[2][4][1]->value) {
    ctr->wf += ctr->locProbs[ctr->fparray[0]-1];
    if (debug) cout << "Adding to wf: p=" << ctr->locProbs[ctr->fparray[0]-1] << " cause=" << ctr->locCauses[ctr->fparray[0]-1] << " loc=" << ctr->fparray[0]-1 << endl;
  }
}


// Obtain weights of the six types of edges. Labeling of edges a, b, c, d, and e is consistent with Fig. 4 in arXiv:1009.3686v1.
void Cube::finalizeRelativeWeights (double p1, double p2, double pm, bool circErr, Counter* ctr) {

  if (circErr) {

    assert (!ctr->usePath && !ctr->countPath);

    if (ctr->wa > 1e-10)
      wa = -1 * log (ctr->wa);
    if (ctr->wb > 1e-10)
      wb = -1 * log (ctr->wb);
    if (ctr->wc > 1e-10)
      wc = -1 * log (ctr->wc);
    if (ctr->wd > 1e-10)
      wd = -1 * log (ctr->wd);
    if (ctr->we > 1e-10)
      we = -1 * log (ctr->we);
    if (ctr->wf > 1e-10)
      wf = -1 * log (ctr->wf);

    if (debug) cout << "getRelativeWeights(p1=" << p1 << " p2=" << p2 << " pm=" << pm << "): " << endl;
    if (debug) cout << "  probs before taking log: wa=" << ctr->wa << " wb=" << ctr->wb << " wc=" << ctr->wc << " wd=" << ctr->wd << " we=" << ctr->we << " wf=" << ctr->wf << endl;
    if (debug) cout << "  calculated initial weights: wa=" << wa << " wb=" << wb << " wc=" << wc << " wd=" << wd << " we=" << we << " wf=" << wf << endl;

    /* Hardcoded values from arXiv:1009.3686v1.
    wa = 1000; // vertical edge along the time dimension
    wb = 1000; // horizontal edge along the x axis
    wc = 1000; // diagonal on a face
    wd = 1000; // horizontal edge along the y axis
    we = 1000; // diagonal on a face
    wf = 1000; // diagonal inside the cube

    // Here we calculate the edge weights according to their probabilities.
    if (p2 > 1e-10 || pm > 1e-10)
      wa = -1 * log (16*p2/15 + pm);
    if (p2 > 1e-10 || p1 > 1e-10)
      wb = -1 * log (8*p2/15 + 4*p1/3);
    if (p2 > 1e-10)
      wc = -1 * log (16*p2/15);
    if (p2 > 1e-10 || p1 > 1e-10)
      wd = -1 * log (32*p2/15 + 4*p1/3);
    if (p2 > 1e-10)
      we = -1 * log (8*p2/15);
    if (p2 > 1e-10)
      wf = -1 * log (8*p2/15);
    */

    vector<double> alt;
    bool change;

    // Here we change the edge weights so that triangle inequality is always satisfied. If an edge is too long and a shorter alternate path exists, the substitution is made.
    do {
      change = false;

      // Substitute for wa
      alt.push_back (wa);
      alt.push_back (wb + wc);
      alt.push_back (wd + we);
      alt.push_back (wb + wd + wf);
      alt.push_back (wc + we + wf);
      if (wa > min_double (alt) + 1e-6) {
        wa = min_double (alt);
        change = true;
      }
      alt.clear ();

      // Substitute for wb
      alt.push_back (wb);
      alt.push_back (wa + wc);
      alt.push_back (wa + wd + wf);
      alt.push_back (wc + wd + we);
      if (wb > min_double (alt) + 1e-6) {
        wb = min_double (alt);
        change = true;
      }

      alt.clear ();

      // Substitute for wc
      alt.push_back (wc);
      alt.push_back (wa + wb);
      alt.push_back (wd + wf);
      alt.push_back (wb + wd + we);
      alt.push_back (wa + we + wf);
      if (wc > min_double (alt) + 1e-6) {
        wc = min_double (alt);
        change = true;
      }
      alt.clear ();

      // Substitute for wd
      alt.push_back (wd);
      alt.push_back (wa + we);
      alt.push_back (wa + wb + wf);
      alt.push_back (wb + wc + we);
      alt.push_back (wc + wf);
      if (wd > min_double (alt) + 1e-6) {
        wd = min_double (alt);
        change = true;
      }
      alt.clear ();

      // Substitute for we
      alt.push_back (we);
      alt.push_back (wa + wd);
      alt.push_back (wb + wc + wd);
      alt.push_back (wa + wc + wf);
      if (we > min_double (alt) + 1e-6) {
        we = min_double (alt);
        change = true;
      }
      alt.clear ();

      // Substitute for wf
      alt.push_back (wf);
      alt.push_back (wa + wb + wd);
      alt.push_back (wc + wd);
      alt.push_back (wa + wc + we);
      if (wf > min_double (alt) + 1e-6) {
        wf = min_double (alt);
        change = true;
      }
      alt.clear ();

    } while (change);
  }

  if (!circErr) {
    if (pm!=0) wa = log ((1 - pm) / pm); else wa = 0; // vertical edge along the time dimension
    if (p1!=0) wb = log ((1 - p1) / p1); else wb = 0; // horizontal edge along the x axis
    if (p1!=0) wd = log ((1 - p1) / p1); else wd = 0; // horizontal edge along the y axis
    wc = wa + wb;      // diagonal on a face
    we = wa + wd;      // diagonal on a face
    wf = wa + wb + wd; // diagonal inside the cube
  }

  if (debug) cout << "  final weights:  " << wa << " " << wb << " " << wc << " " << wd << " " << we << " " << wf << endl;

}


// For each location in the cube, initialize the pointers indicating the nearest vertex (non-tirvial syndrome) in five directions in the cube. This structure is used for prunning unnecesary edges when creating graph for matching.
void Cube::initNearest (void){

  // Populate nearestUp
  for (int t = 0; t < TSize; t++) {
    for (int i = 0; i < xSize; i++) {
      Vertex* pPrev = NULL;
      for (int j = 0; j < ySize; j++) {
        if (! ((i%2 && j%2) || (!(i%2) && !(j%2))) )
          continue;
        if (syndromes[i][j][t]->value)
          pPrev = syndromes[i][j][t]->pVert;
      }     
      for (int j = 0; j < ySize; j++) {
        if (! ((i%2 && j%2) || (!(i%2) && !(j%2))) )
          continue;
        syndromes[i][j][t]->nearestUp = pPrev;
        if (syndromes[i][j][t]->value)
          pPrev = syndromes[i][j][t]->pVert;
      }
    }
  }

  // Populate nearestDown
  for (int t = 0; t < TSize; t++) {
    for (int i = 0; i < xSize; i++) {
      Vertex* pPrev = NULL;
      for (int j = ySize - 1; j >= 0; j--) {
        if (! ((i%2 && j%2) || (!(i%2) && !(j%2))) )
          continue;
        if (syndromes[i][j][t]->value)
          pPrev = syndromes[i][j][t]->pVert;
      }
      for (int j = ySize - 1; j >= 0; j--) {
        if (! ((i%2 && j%2) || (!(i%2) && !(j%2))) )
          continue;
        syndromes[i][j][t]->nearestDown = pPrev;
        if (syndromes[i][j][t]->value)
          pPrev = syndromes[i][j][t]->pVert;
      }
    }
  }

  // Populate nearestLeft
  for (int t = 0; t < TSize; t++) {
    for (int j = 0; j < ySize; j++) {
      Vertex* pPrev = NULL;
      for (int i = 0; i < xSize; i++) {
        if (! ((i%2 && j%2) || (!(i%2) && !(j%2))) )
          continue;
        if (syndromes[i][j][t]->value)
          pPrev = syndromes[i][j][t]->pVert;
      }
      for (int i = 0; i < xSize; i++) {
        if (! ((i%2 && j%2) || (!(i%2) && !(j%2))) )
          continue;
        syndromes[i][j][t]->nearestLeft = pPrev;
        if (syndromes[i][j][t]->value)
          pPrev = syndromes[i][j][t]->pVert;
      }
    }
  }

  // Populate nearestRight
  for (int t = 0; t < TSize; t++) {
    for (int j = 0; j < ySize; j++) {
      Vertex* pPrev = NULL;
      for (int i = xSize - 1; i >= 0; i--) {
        if (! ((i%2 && j%2) || (!(i%2) && !(j%2))) )
          continue;
        if (syndromes[i][j][t]->value)
          pPrev = syndromes[i][j][t]->pVert;
      }
      for (int i = xSize - 1; i >= 0; i--) {
        if (! ((i%2 && j%2) || (!(i%2) && !(j%2))) )
          continue;
        syndromes[i][j][t]->nearestRight = pPrev;
        if (syndromes[i][j][t]->value)
          pPrev = syndromes[i][j][t]->pVert;
      }
    }
  }

  // Populate nearestAfter
  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      if (! ((i%2 && j%2) || (!(i%2) && !(j%2))) )
        continue;
      Vertex* pPrev = NULL;
      for (int t = TSize - 1; t >= 0; t--) {
        if (syndromes[i][j][t]->value)
          pPrev = syndromes[i][j][t]->pVert;
      }
      for (int t = TSize - 1; t >= 0; t--) {
        syndromes[i][j][t]->nearestAfter = pPrev;
        if (syndromes[i][j][t]->value)
          pPrev = syndromes[i][j][t]->pVert;
      }
    }
  }

  // Debugging output
  if (debug) {
    for (int i = 0; i < xSize; i++) {
      for (int j = 0; j < ySize; j++) {
        if (! ((i%2 && j%2) || (!(i%2) && !(j%2))) )
          continue;
        for (int t = 0; t < TSize; t++) {
          cout << "nearest[" << i << "][" << j << "][" << t << "]:";
          if (syndromes[i][j][t]->nearestUp != NULL) cout << " Up=" << syndromes[i][j][t]->nearestUp->y;
          if (syndromes[i][j][t]->nearestDown != NULL) cout << " Down=" << syndromes[i][j][t]->nearestDown->y;
          if (syndromes[i][j][t]->nearestRight != NULL) cout << " Right=" << syndromes[i][j][t]->nearestRight->x;
          if (syndromes[i][j][t]->nearestLeft != NULL) cout << " Left=" << syndromes[i][j][t]->nearestLeft->x;
          if (syndromes[i][j][t]->nearestAfter != NULL) cout << " After=" << syndromes[i][j][t]->nearestAfter->t;
          cout << endl; 
        }
      }
    }
  }

}


// Obtain candidate vertices for matching using the complete XXXX and ZZZZ syndrome history.
void Cube::getVertices (vector<Vertex*> & verticesXXXX, vector<Vertex*> & verticesZZZZ){

  int xid = 0, zid = 0;

  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      if (! ((i%2 && j%2) || (!(i%2) && !(j%2))) )
        continue;
      bool prev = false;
      for (int t = 0; t < TSize; t++) {
        bool newVert = false;
        if (syndromes[i][j][t]->value != prev && !lkgHistory[i][j][t])
          newVert = true;
        if (!lkgHistory[i][j][t])
          prev = syndromes[i][j][t]->value;
        if (!newVert)
          syndromes[i][j][t]->value = false;
        if (newVert) {
          syndromes[i][j][t]->value = true;
          if (i%2 && j%2) {
            Vertex* pVertex = new Vertex (i, j, t, zid++, syndromes[i][j][t]->syndromeID);
            verticesZZZZ.push_back(pVertex);
            syndromes[i][j][t]->pVert = pVertex;
            if (debug) cout << "Created new Z vertex: i=" << i << " j=" << j << " t=" << t << endl;
          } else {
            Vertex* pVertex = new Vertex (i, j, t, xid++, syndromes[i][j][t]->syndromeID);
            verticesXXXX.push_back(pVertex);
            syndromes[i][j][t]->pVert = pVertex;
            if (debug) cout << "Created new X vertex: i=" << i << " j=" << j << " t=" << t << endl;
          }
        }
      }
    }
  }

}


// Conncet the vertices that need to be matched by edges. The result in stored in the vector edges.
void Cube::getEdges (vector<Vertex*> & vertices, vector<MyEdge*> & edges){

  // For each vertex u, the body of the loop finds all other vertices v that should be connected by an edge. Edge prunning inspired by arXiv 1307:1740 is used -- a sufficiently large rectangle is formed around vertex v and then we connect to all vertices inside this rectangle.
  for (int u = 0; u < (int)vertices.size(); u++) {
    int ux = vertices[u]->x;
    int uy = vertices[u]->y;
    int ut = vertices[u]->t;
    if (debug) cout << "New vertex u: ux=" << ux << " uy="<< uy << " ut=" << ut << endl;

    assert (syndromes[ux][uy][ut]->nearestLeft != NULL);
    assert (syndromes[ux][uy][ut]->nearestUp != NULL);

    int vLeft;
    int vRight;
    int vUp;
    int vDown;
    
    vLeft = syndromes[ux][uy][ut]->nearestLeft->x;
    vRight = syndromes[ux][uy][ut]->nearestRight->x;

    if (vLeft != vRight) {
      vRight = (vRight + 2) % xSize;
    }  
    
    vUp = syndromes[ux][uy][ut]->nearestUp->y;
    vDown = syndromes[ux][uy][ut]->nearestDown->y;

    if (vUp != vDown) {
      vDown = (vDown + 2) % ySize;
    }

    int i = vLeft; do {
      int j = vUp; do {
        Vertex* pv;
        if (syndromes[i%xSize][j%ySize][ut]->value && syndromes[i%xSize][j%ySize][ut]->pVert->id != u)
          pv = syndromes[i%xSize][j%ySize][ut]->pVert;
        else
          pv = syndromes[i%xSize][j%ySize][ut]->nearestAfter;

        if (pv == NULL || pv->id == u) {
          j += 2;
          continue;
        }
        int vx  = pv->x;
        int vy  = pv->y;
        int vt  = pv->t;
        int vid = pv->id;

        // Calculate the appropriate distance of the vertex pair and add the edge
        double dist = getDistance (ux, uy, ut, vx, vy, vt);
        MyEdge* newEdge = new MyEdge (u,vid,dist);
        edges.push_back(newEdge);
        if (debug) cout << "Adding new edge: u=" << u << "(" << ux << "," << uy << "," << ut << ") and vid=" << vid << "(" << vx << "," << vy << "," << vt << "): dist=" << dist <<  endl;
        j += 2;
      } while ((j%ySize) != vDown);
      i+= 2;
    } while ((i%xSize) != vRight);
  }
}


// Conncet the vertices that need to be matched by edges. The result in stored in the vector edges.
void Cube::getEdgesLkgDetBoost (vector<Vertex*> & vertices, vector<MyEdge*> & edges, pauli errType){

  /*////////// BOOST EXAMPLE //////////
  // Create undirected weighted graph
  typedef adjacency_list<listS, vecS, undirectedS, no_property, property<edge_weight_t, double> > mygraph;
  mygraph g;
  add_edge(0, 1, 0.5, g);
  add_edge(1, 3, 2.0, g);
  add_edge(1, 4, 0.1, g);
  add_edge(0, 2, 2.5, g);
  add_edge(2, 4, 0.2, g);
  add_edge(2, 5, 3.0, g);
  mygraph::vertex_descriptor vs = vertex(0, g);

  // Run Dijkstra
  vector<mygraph::vertex_descriptor> parents(num_vertices(g)); // To store parents
  vector<double> distances(num_vertices(g)); // To store distances
  dijkstra_shortest_paths(g, vs, predecessor_map(&parents[0]).distance_map(&distances[0]));

  // Report distance 0->5
  mygraph::vertex_descriptor vd = vertex(5, g);
  cout << "distance(5) = " << distances[5] << ", ";

  // Report all results
  cout << "distances and parents:" << endl;
  graph_traits < mygraph >::vertex_iterator vertexIterator, vend;
  for (boost::tie(vertexIterator, vend) = boost::vertices(g); vertexIterator != vend; ++vertexIterator) 
    {
      cout << "distance(" << *vertexIterator << ") = " << distances[*vertexIterator] << ", ";
      cout << "parent(" << *vertexIterator << ") = " << parents[*vertexIterator] << endl;
    }
  cout << endl;
  ////////// END OF BOOST EXAMPLE //////////*/

  // Prepare empty undirected weighted graph
  mygraph g(numSynd);


  // Initialize edge distances
  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      if ((i%2 && !(j%2)) || (!(i%2) && j%2))
        continue;
      if (errType == X && !(i%2) && !(j%2))
        continue;
      if (errType == Z && i%2 && j%2)
        continue;
      for (int t = 0; t < TSize; t++) {
        int this_id = syndromes[i][j][t]->syndromeID;
        // a edge
        if (t<TSize-1) {
          int top_id = syndromes[i][j][t+1]->syndromeID;
          // The weight wa should be 0 in case of a leakage
          double thisWa = wa;
          if (lkgHistory[i][j][t])
            thisWa = 0;
          updateEdge(this_id, top_id, thisWa, g);
          if (debug) cout << "New edge (a)" << i << "," << j << "," << t << "<->" << i << "," << j << "," << t+1 << " wt. " << thisWa << endl;
        }
        // b edge
        int north_id = syndromes[i][(j-2+ySize)%ySize][t]->syndromeID;
        updateEdge(this_id, north_id, wb, g);
        if (debug) cout << "New edge (b)" << i << "," << j << "," << t << "<->" << i << "," << (j-2+ySize)%ySize << "," << t << " wt. " << wb << endl;
        // c edge
        if (t<TSize-1) {
          int southtop_id = syndromes[i][(j+2)%ySize][t+1]->syndromeID;
          updateEdge(this_id, southtop_id, wc, g);
          if (debug) cout << "New edge (c)" << i << "," << j << "," << t << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << wc << endl;
        }
        // d edge
        int east_id = syndromes[(i+2)%xSize][j][t]->syndromeID;
        updateEdge(this_id, east_id, wd, g);
        if (debug) cout << "New edge (d)" << i << "," << j << "," << t << "<->" << (i+2)%xSize << "," << j << "," << t << " wt. " << wd << endl;
        // e edge
        if (t<TSize-1) {
          int easttop_id = syndromes[(i+2)%xSize][j][t+1]->syndromeID;
          updateEdge(this_id, easttop_id, we, g);
          if (debug) cout << "New edge (e)" << i << "," << j << "," << t << "<->" << (i+2)%xSize << "," << j << "," << t+1 << " wt. " << we << endl;
        }
        // f edge
        if (t<TSize-1) {
          int swtop_id = syndromes[(i-2+xSize)%xSize][(j+2)%ySize][t+1]->syndromeID;
          updateEdge(this_id, swtop_id, wf, g);
          if (debug) cout << "New edge (f)" << i << "," << j << "," << t << "<->" << (i-2+xSize)%xSize << "," << (j+2)%ySize << "," << t+1 << " wt. " << wf << endl;
        }
      }
    }
  }

  // Update edge disgtances for current leakage history - CIRCUIT MODEL
  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      for (int t = 0; t < TSize; t++) {
        if (!lkgHistory[i][j][t])
          continue;
        if (debug && corLeak == circuit) cout << "Updating edges for leakage at i,j,t: " << i << " " << j << " " << t << endl;	
        // Add edges - circuit model, ancilla
	if ((t <= TSize - 2 && corLeak == circuit && errType == Z && (!(i%2) && !(j%2))) || 
	    (t <= TSize - 2 && corLeak == circuit && errType == X && ((i%2) && (j%2)))) {
          int tc_alpha_id = syndromes[i][(j-2+ySize)%ySize][t]->syndromeID;
          int tc_beta_id  = syndromes[(i-2+xSize)%xSize][j][t]->syndromeID;
          int tc_star_id  = syndromes[i][j][t]->syndromeID;
          int tn_star_id  = syndromes[i][j][t+1]->syndromeID;
          int tn_gamma_id = syndromes[(i+2)%xSize][j][t+1]->syndromeID;
          int tn_delta_id = syndromes[i][(j+2)%ySize][t+1]->syndromeID;
          // Edge w. p. 1/2*1/5
          updateEdge(tc_star_id, tc_alpha_id, -1*log((double)1/10), g);
          if (debug) cout << "New edge " << i << "," << j << "," << t+1 << "<->" << i << "," << (j-2+ySize)%ySize << "," << t << " wt. " << -1*log((double)1/10) << endl;
          // Edge w. p. 1/2*2/5
          updateEdge(tc_star_id, tc_beta_id, -1*log((double)2/10), g);
          if (debug) cout << "New edge " << i << "," << j << "," << t+1 << "<->" << (i-2+xSize)%xSize << "," << j << "," << t << " wt. " << -1*log((double)2/10) << endl;
          // Edge w. p. 1/2*3/5
          updateEdge(tn_star_id, tn_gamma_id, -1*log((double)3/10), g);
          if (debug) cout << "New edge " << i << "," << j << "," << t+1 << "<->" << (i+2)%xSize << "," << j << "," << t+1 << " wt. " << -1*log((double)3/10) << endl;
          // Edge w. p. 1/2*4/5
          updateEdge(tn_star_id, tn_delta_id, -1*log((double)4/10), g);
          if (debug) cout << "New edge " << i << "," << j << "," << t+1 << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)4/10) << endl;
        }

        // Add edges - circuit model, ancilla
        if ((corLeak == circuit && errType == Z && ((i%2) && (j%2))) ||
            (corLeak == circuit && errType == X && (!(i%2) && !(j%2)))) {
          int tc_1_id  = syndromes[(i-1+xSize)%xSize][(j-1+ySize)%ySize][t]->syndromeID;
          int tc_2_id  = syndromes[(i+1)%xSize][(j-1+ySize)%ySize][t]->syndromeID;
          // Edge w. p. 1/2*1/5
          updateEdge(tc_1_id, tc_2_id, -1*log((double)1/10), g);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j-1+ySize)%ySize << "," << t << "<->" << (i+1)%xSize << "," << (j-1+ySize)%ySize << "," << t << " wt. " << -1*log((double)1/10) << endl;
          if (t >= TSize - 1)
            break;
          int tn_3_id = syndromes[(i-1+xSize)%xSize][(j+1)%ySize][t+1]->syndromeID;
          int tn_4_id = syndromes[(i+1)%xSize][(j+1)%ySize][t+1]->syndromeID;
          // Edge w. p. 1/2*2/5
          updateEdge(tc_1_id, tn_3_id, -1*log((double)2/10), g);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j-1+ySize)%ySize << "," << t << "<->" << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)2/10) << endl;
          // Edge w. p. 1/2*3/5
          updateEdge(tc_2_id, tn_4_id, -1*log((double)3/10), g);
          if (debug) cout << "New edge " << (i+1)%xSize << "," << (j-1+ySize)%ySize << "," << t << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)3/10) << endl;
          // Edge w. p. 1/2*4/5
          updateEdge(tn_3_id, tn_4_id, -1*log((double)4/10), g);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t+1 << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)4/10) << endl;
        }
        // Add edges - circuit model, data
        if ((t <= TSize - 2 && corLeak == circuit && errType == Z && (!(i%2) && (j%2))) ||
            (t <= TSize - 2 && corLeak == circuit && errType == X && ((i%2) && !(j%2)))) {
          int tc_alpha_id = syndromes[i][(j-1+ySize)%ySize][t]->syndromeID;
          int tc_beta_id  = syndromes[i][(j+1)%ySize][t]->syndromeID;
          int tc_east_id  = syndromes[(i+2)%xSize][(j-1+ySize)%ySize][t]->syndromeID;
          int tn_alpha_id = syndromes[i][(j-1+ySize)%ySize][t+1]->syndromeID;
          int tn_beta_id  = syndromes[i][(j+1)%ySize][t+1]->syndromeID;
          int tn_west_id  = syndromes[(i-2+xSize)%xSize][(j+1)%ySize][t+1]->syndromeID;

          // Edge w. p. 1/2*3/9
          updateEdge(tc_beta_id, tn_beta_id, -1*log((double)3/18), g);
          if (debug) cout << "New edge " << i << "," << (j+1)%ySize << "," << t+1 << "<->" << i << "," << (j+1)%ySize << "," << t << " wt. " << -1*log((double)3/18) << endl;
          // Edge w. p. 1/2*4/9
          updateEdge(tc_east_id, tn_beta_id, -1*log((double)4/18), g);
          if (debug) cout << "New edge " << i << "," << (j+1)%ySize << "," << t+1 << "<->" << (i+2)%xSize << "," << (j-1+ySize)%ySize << "," << t << " wt. " << -1*log((double)4/18) << endl;
          // Edge w. p. 1/2*5/9
          updateEdge(tn_beta_id, tn_west_id, -1*log((double)5/18), g);
          if (debug) cout << "New edge " << i << "," << (j+1)%ySize << "," << t+1 << "<->" << (i-2+xSize)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)5/18) << endl;
          // Edge w. p. 1/2*6/9
	  updateEdge(tc_alpha_id, tn_alpha_id, -1*log((double)6/18), g);
          if (debug) cout << "New edge " << i << "," << (j-1+ySize)%ySize << "," << t << "<->" << i << "," << (j-1+ySize)%ySize << "," << t+1 << " wt. " << -1*log((double)6/18) << endl;
          // Edge w. p. 1/2*8/9
          updateEdge(tn_beta_id, tn_alpha_id, -1*log((double)8/18), g);
          if (debug) cout << "New edge " << i << "," << (j+1)%ySize << "," << t+1 << "<->" << i << "," << (j-1+ySize)%ySize << "," << t+1 << " wt. " << -1*log((double)8/18) << endl;
        }

        // Add edges - circuit model, data
        if ((corLeak == circuit && errType == Z && ((i%2) && !(j%2))) ||
            (corLeak == circuit && errType == X && (!(i%2) && (j%2)))) {
          int tc_1_id  = syndromes[(i-1+xSize)%xSize][j][t]->syndromeID;
          int tc_2_id  = syndromes[(i+1)%xSize][j][t]->syndromeID;
          // Edge w. p. 1/2*3/9
          updateEdge(tc_1_id, tc_2_id, -1*log((double)3/18), g);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << j << "," << t << "<->" << (i+1)%xSize << "," << j << "," << t << " wt. " << -1*log((double)3/18) << endl;
          if (t >= TSize - 1)
            break;
          int tn_1_id = syndromes[(i-1+xSize)%xSize][j][t+1]->syndromeID;
          int tn_2_id = syndromes[(i+1)%xSize][j][t+1]->syndromeID;
          // Edge w. p. 1/2*4/9
          updateEdge(tc_2_id, tn_2_id, -1*log((double)4/18), g);
          if (debug) cout << "New edge " << (i+1)%xSize << "," << j << "," << t << "<->" << (i+1)%xSize << "," << j << "," << t+1 << " wt. " << -1*log((double)4/18) << endl;
          // Edge w. p. 1/2*5/9
          updateEdge(tc_1_id, tn_1_id, -1*log((double)5/18), g);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << j << "," << t << "<->" << (i-1+xSize)%xSize << "," << j << "," << t+1 << " wt. " << -1*log((double)5/18) << endl;
          // Edge w. p. 1/2*8/9
          updateEdge(tn_1_id, tn_2_id, -1*log((double)8/18), g);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << j << "," << t+1 << "<->" << (i+1)%xSize << "," << j << "," << t+1 << " wt. " << -1*log((double)8/18) << endl;
        }
      }
    }
  }


  // Update edge disgtances for current leakage history - QUICK MODEL
  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      for (int t = 0; t < TSize; t++) {
        if ((i%2 && !(j%2)) || (!(i%2) && j%2) || !lkgHistory[i][j][t])
          continue;
        if (debug && corLeak == quick) cout << "Updating edges for leakage at i,j,t: " << i << " " << j << " " << t << endl;
        // Add edges - quick model
	if ((corLeak == quick && errType == Z && ((i%2) && (j%2))) || 
	    (corLeak == quick && errType == X && (!(i%2) && !(j%2)))) {
          int tc_3_id = syndromes[(i-1+xSize)%xSize][(j+1)%ySize][t]->syndromeID;
          int tc_4_id = syndromes[(i+1)%xSize][(j+1)%ySize][t]->syndromeID;
          // Edge w. p. 1/2*6/11
          updateEdge(tc_3_id, tc_4_id, -1*log((double)6/22), g);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t << " wt. " << -1*log((double)6/22) << endl;
          if (t >= 1) {
            int tp_1_id  = syndromes[(i-1+xSize)%xSize][(j-1+ySize)%ySize][t-1]->syndromeID;
            int tp_2_id  = syndromes[(i+1)%xSize][(j-1+ySize)%ySize][t-1]->syndromeID;
            // Edge w. p. 1/2*1/11
            updateEdge(tp_1_id, tp_2_id, -1*log((double)1/22), g);
            if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j-1+ySize)%ySize << "," << t-1 << "<->" << (i+1)%xSize << "," << (j-1+ySize)%ySize << "," << t-1 << " wt. " << -1*log((double)1/22) << endl;
            // Edge w. p. 1/2*2/11
            updateEdge(tp_1_id, tc_3_id, -1*log((double)2/22), g);
            if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j-1+ySize)%ySize << "," << t-1 << "<->" << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t << " wt. " << -1*log((double)2/22) << endl;
            // Edge w. p. 1/2*3/11
            updateEdge(tp_2_id, tc_4_id, -1*log((double)3/22), g);
            if (debug) cout << "New edge " << (i+1)%xSize << "," << (j-1+ySize)%ySize << "," << t-1 << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t << " wt. " << -1*log((double)3/22) << endl;
	  }
          if (t <= TSize - 2) {
            int tn_3_id = syndromes[(i-1+xSize)%xSize][(j+1)%ySize][t+1]->syndromeID;
            int tn_4_id = syndromes[(i+1)%xSize][(j+1)%ySize][t+1]->syndromeID;
            // Edge w. p. 1/2*7/11
	    updateEdge(tc_4_id, tn_4_id, -1*log((double)7/22), g);
            if (debug) cout << "New edge " << (i+1)%xSize << "," << (j+1)%ySize << "," << t << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)7/22) << endl;
            // Edge w. p. 1/2*8/11
            updateEdge(tc_3_id, tn_3_id, -1*log((double)8/22), g);
            if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t << "<->" << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)8/22) << endl;
            // Edge w. p. 1/2*10/11
            updateEdge(tn_3_id, tn_4_id, -1*log((double)10/22), g);
            if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t+1 << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)10/22) << endl;
          }
        }

        // Add edges - quick model
	if ((corLeak == quick && errType == Z && (!(i%2) && !(j%2))) ||
            (corLeak == quick && errType == X && ((i%2) && (j%2)))) {
          int tp_alpha_id = 0, tp_beta_id = 0, tp_star_id = 0;
          if (t >= 1) {
            tp_alpha_id       = syndromes[i][(j-2+ySize)%ySize][t-1]->syndromeID;
            tp_beta_id        = syndromes[(i-2+xSize)%xSize][j][t-1]->syndromeID;
            tp_star_id        = syndromes[i][j][t-1]->syndromeID;
          }
          int tc_star_id      = syndromes[i][j][t]->syndromeID;
          int tc_gamma_id     = syndromes[(i+2)%xSize][j][t]->syndromeID;
          int tc_delta_id     = syndromes[i][(j+2)%ySize][t]->syndromeID;
          int tn_delta_id = 0, tn_star_id = 0, tn_deltawest_id = 0;
          if (t <= TSize - 2) {
            tn_delta_id      = syndromes[i][(j+2)%ySize][t+1]->syndromeID;
            tn_star_id       = syndromes[i][j][t+1]->syndromeID;
            tn_deltawest_id  = syndromes[(i-2+xSize)%xSize][(j+2)%ySize][t+1]->syndromeID;
          }
          if (t >= 1 && t <= TSize - 2) {
            // Edge w. p. 1/2*1/11
            updateEdge(tp_alpha_id, tn_delta_id, -1*log((double)1/22), g);
            if (debug) cout << "New edge " << i << "," << (j-2+ySize)%ySize << "," << t-1 << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)1/22) << endl;
            // Edge w. p. 1/2*2/11
            updateEdge(tp_beta_id, tn_delta_id, -1*log((double)2/22), g);
            if (debug) cout << "New edge " << (i-2+xSize)%xSize << "," << j << "," << t-1 << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)2/22) << endl;
          }
          if (t >= 1 && t <= TSize - 1) {
            // Edge w. p. 1/2*5/11
            updateEdge(tp_star_id, tc_star_id, -1*log((double)5/22), g);
            if (debug) cout << "New edge " << i << "," << j << "," << t-1 << "<->" << i << "," << j << "," << t+1 << " wt. " << -1*log((double)5/22) << endl;
          }
          if (t <= TSize - 2) {
            // Edge w. p. 1/2*6/11
            updateEdge(tc_delta_id, tn_delta_id, -1*log((double)6/22), g);
            if (debug) cout << "New edge " << i << "," << (j+2)%ySize << "," << t << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)6/22) << endl;
            // Edge w. p. 1/2*(7+3)/11 - 3/22*7/22
            updateEdge(tc_gamma_id, tn_delta_id, -1*log((double)199/484), g);
            if (debug) cout << "New edge " << (i+2)%xSize << "," << j << "," << t << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)199/484) << endl;
            // Edge w. p. 1/2*8/11
            updateEdge(tn_delta_id, tn_deltawest_id, -1*log((double)8/22), g);
            if (debug) cout << "New edge " << i << "," << (j+2)%ySize << "," << t+1 << "<->" << (i-2+xSize)%xSize << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)8/22) << endl;
            // Edge w. p. 1/2*10/11
            updateEdge(tn_delta_id, tc_star_id, -1*log((double)10/22), g);
            if (debug) cout << "New edge " << i << "," << (j+2)%ySize << "," << t+1 << "<->" << i << "," << j << "," << t+1 << " wt. " << -1*log((double)10/22) << endl;
          }
        }
      }
    }
  }

  // For each vertex u, connect to every other vertex v by an edge. No edge prunning performed.
  for (int u = 0; u < (int)vertices.size(); u++) {
    int uid = vertices[u]->id;
    int usid = vertices[u]->syndromeID;
    // Run Dijkstra from vertex u
    mygraph::vertex_descriptor vs = vertex(usid, g);
    vector<mygraph::vertex_descriptor> parents(num_vertices(g)); // To store parents
    vector<double> distances(num_vertices(g)); // To store distances
    dijkstra_shortest_paths(g, vs, predecessor_map(&parents[0]).distance_map(&distances[0]));
    for (int v = 0; v < (int)vertices.size(); v++) {
      if (u == v)
        continue;
      // Add the edge with its distance
      int vid = vertices[v]->id;
      int vsid = vertices[v]->syndromeID;
      double dist = distances[vsid];
      MyEdge* newEdge = new MyEdge (uid,vid,dist);
      edges.push_back(newEdge);
      if (debug) {
        int ux  = vertices[u]->x;
        int uy  = vertices[u]->y;
        int ut  = vertices[u]->t;
        int vx  = vertices[v]->x;
        int vy  = vertices[v]->y;
        int vt  = vertices[v]->t;
        cout << "Adding new edge: uid=" << uid << "(" << ux << "," << uy << "," << ut << ") and vid=" << vid << "(" << vx << "," << vy << "," << vt << "): dist=" << dist <<  endl;
      }
    }
  }
}


// Find out if an edge in the Boost library exists, and update its weight.
void Cube::updateEdge(int v1, int v2, double wtAdd, mygraph & g) {

  // Find out if the edge exists and read it into ed
  pair<graph_traits<mygraph>::edge_descriptor, bool> ed = edge(v1, v2, g);

  if (!ed.second) {
    add_edge(v1 , v2, wtAdd, g);
  } else {
    double wtOld = get(edge_weight_t(), g, ed.first);
    double wtA = exp(-1.0*wtAdd);
    double wtB = exp(-1.0*wtOld);
    put(edge_weight_t(), g, ed.first, -1.0*log(wtA + wtB -wtA*wtB));
  }

}


// Conncet the vertices that need to be matched by edges. The result in stored in the vector edges.
void Cube::getEdgesLkgDetFloydWarshall (vector<Vertex*> & vertices, vector<MyEdge*> & edges, pauli errType){

  // Implements Floyd Warshall algorithm
  double MDist[numSynd][numSynd];
  for (int i = 0; i < numSynd; i++) {
    for (int j = 0; j < numSynd; j++) {
      MDist[i][j] = INT_MAX;
    }
  }

  // Initialize edge distances
  for (int i = 0; i < numSynd; i++)
    MDist[i][i] = 0;

  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      if ((i%2 && !(j%2)) || (!(i%2) && j%2))
        continue;
      for (int t = 0; t < TSize; t++) {
        int this_id = syndromes[i][j][t]->syndromeID;
        // a edge
        if (t<TSize-1) {
          int top_id = syndromes[i][j][t+1]->syndromeID;
          MDist[this_id][top_id] = wa;
          MDist[top_id][this_id] = wa;
          if (debug) cout << "New edge (a)" << i << "," << j << "," << t << "<->" << i << "," << j << "," << t+1 << " wt. " << wa << endl;
        }
        // b edge
        int north_id = syndromes[i][(j-2+ySize)%ySize][t]->syndromeID;
        MDist[this_id][north_id] = wb;
        MDist[north_id][this_id] = wb;
        if (debug) cout << "New edge (b)" << i << "," << j << "," << t << "<->" << i << "," << (j-2+ySize)%ySize << "," << t << " wt. " << wb << endl;
        // c edge
        if (t<TSize-1) {
          int southtop_id = syndromes[i][(j+2)%ySize][t+1]->syndromeID;
          MDist[this_id][southtop_id] = wc;
          MDist[southtop_id][this_id] = wc;
          if (debug) cout << "New edge (c)" << i << "," << j << "," << t << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << wc << endl;
        }
        // d edge
        int east_id = syndromes[(i+2)%xSize][j][t]->syndromeID;
        MDist[this_id][east_id] = wd;
        MDist[east_id][this_id] = wd;
        if (debug) cout << "New edge (d)" << i << "," << j << "," << t << "<->" << (i+2)%xSize << "," << j << "," << t << " wt. " << wd << endl;
        // e edge
        if (t<TSize-1) {
          int easttop_id = syndromes[(i+2)%xSize][j][t+1]->syndromeID;
          MDist[this_id][easttop_id] = we;
          MDist[easttop_id][this_id] = we;
          if (debug) cout << "New edge (e)" << i << "," << j << "," << t << "<->" << (i+2)%xSize << "," << j << "," << t+1 << " wt. " << we << endl;
        }
        // f edge
        if (t<TSize-1) {
          int swtop_id = syndromes[(i-2+xSize)%xSize][(j+2)%ySize][t+1]->syndromeID;
          MDist[this_id][swtop_id] = wf;
          MDist[swtop_id][this_id] = wf;
          if (debug) cout << "New edge (f)" << i << "," << j << "," << t << "<->" << (i-2+xSize)%xSize << "," << (j+2)%ySize << "," << t+1 << " wt. " << wf << endl;
        }
      }
    }
  }

  // Update edge disgtances for current leakage history - CIRCUIT MODEL
  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      for (int t = 0; t < TSize; t++) {
        if (!lkgHistory[i][j][t])
          continue;
        if (debug && corLeak == circuit) cout << "Updating edges for leakage at i,j,t: " << i << " " << j << " " << t << endl;	
        // Add edges - circuit model, ancilla
	if ((t <= TSize - 2 && corLeak == circuit && errType == Z && (!(i%2) && !(j%2))) || 
	    (t <= TSize - 2 && corLeak == circuit && errType == X && ((i%2) && (j%2)))) {
          int tc_alpha_id = syndromes[i][(j-2+ySize)%ySize][t]->syndromeID;
          int tc_beta_id  = syndromes[(i-2+xSize)%xSize][j][t]->syndromeID;
          int tc_star_id  = syndromes[i][j][t]->syndromeID;
          int tn_star_id  = syndromes[i][j][t+1]->syndromeID;
          int tn_gamma_id = syndromes[(i+2)%xSize][j][t+1]->syndromeID;
          int tn_delta_id = syndromes[i][(j+2)%ySize][t+1]->syndromeID;
          // Edge w. p. 1/2*1/5
          MDist[tc_star_id][tc_alpha_id] = -1*log((double)1/10);
          MDist[tc_alpha_id][tc_star_id] = -1*log((double)1/10);
          if (debug) cout << "New edge " << i << "," << j << "," << t+1 << "<->" << i << "," << (j-2+ySize)%ySize << "," << t << " wt. " << -1*log((double)1/10) << endl;
          // Edge w. p. 1/2*2/5
          MDist[tc_star_id][tc_beta_id] = -1*log((double)2/10);
          MDist[tc_beta_id][tc_star_id] = -1*log((double)2/10);
          if (debug) cout << "New edge " << i << "," << j << "," << t+1 << "<->" << (i-2+xSize)%xSize << "," << j << "," << t << " wt. " << -1*log((double)2/10) << endl;
          // Edge w. p. 1/2*3/5
          MDist[tn_star_id][tn_gamma_id] = -1*log((double)3/10);
          MDist[tn_gamma_id][tn_star_id] = -1*log((double)3/10);
	  if (debug) cout << "New edge " << i << "," << j << "," << t+1 << "<->" << (i+2)%xSize << "," << j << "," << t+1 << " wt. " << -1*log((double)3/10) << endl;
          // Edge w. p. 1/2*4/5
          MDist[tn_star_id][tn_delta_id] = -1*log((double)4/10);
          MDist[tn_delta_id][tn_star_id] = -1*log((double)4/10);
          if (debug) cout << "New edge " << i << "," << j << "," << t+1 << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)4/10) << endl;
        }

        // Add edges - circuit model, ancilla
        if ((corLeak == circuit && errType == Z && ((i%2) && (j%2))) ||
            (corLeak == circuit && errType == X && (!(i%2) && !(j%2)))) {
          int tc_1_id  = syndromes[(i-1+xSize)%xSize][(j-1+ySize)%ySize][t]->syndromeID;
          int tc_2_id  = syndromes[(i+1)%xSize][(j-1+ySize)%ySize][t]->syndromeID;
          // Edge w. p. 1/2*1/5
          MDist[tc_1_id][tc_2_id] = -1*log((double)1/10);
          MDist[tc_2_id][tc_1_id] = -1*log((double)1/10);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j-1+ySize)%ySize << "," << t << "<->" << (i+1)%xSize << "," << (j-1+ySize)%ySize << "," << t << " wt. " << -1*log((double)1/10) << endl;
          if (t >= TSize - 1)
            break;
          int tn_3_id = syndromes[(i-1+xSize)%xSize][(j+1)%ySize][t+1]->syndromeID;
          int tn_4_id = syndromes[(i+1)%xSize][(j+1)%ySize][t+1]->syndromeID;
          // Edge w. p. 1/2*2/5
          MDist[tc_1_id][tn_3_id] = -1*log((double)2/10);
          MDist[tn_3_id][tc_1_id] = -1*log((double)2/10);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j-1+ySize)%ySize << "," << t << "<->" << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)2/10) << endl;
          // Edge w. p. 1/2*3/5
          MDist[tc_2_id][tn_4_id] = -1*log((double)3/10);
          MDist[tn_4_id][tc_2_id] = -1*log((double)3/10);
          if (debug) cout << "New edge " << (i+1)%xSize << "," << (j-1+ySize)%ySize << "," << t << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)3/10) << endl;
          // Edge w. p. 1/2*4/5
          MDist[tn_3_id][tn_4_id] = -1*log((double)4/10);
          MDist[tn_4_id][tn_3_id] = -1*log((double)4/10);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t+1 << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)4/10) << endl;
        }
        // Add edges - circuit model, data
        if ((t <= TSize - 2 && corLeak == circuit && errType == Z && (!(i%2) && (j%2))) ||
            (t <= TSize - 2 && corLeak == circuit && errType == X && ((i%2) && !(j%2)))) {
          int tc_alpha_id = syndromes[i][(j-1+ySize)%ySize][t]->syndromeID;
          int tc_beta_id  = syndromes[i][(j+1)%ySize][t]->syndromeID;
          int tc_east_id  = syndromes[(i+2)%xSize][(j-1+ySize)%ySize][t]->syndromeID;
          int tn_alpha_id = syndromes[i][(j-1+ySize)%ySize][t+1]->syndromeID;
          int tn_beta_id  = syndromes[i][(j+1)%ySize][t+1]->syndromeID;
          int tn_west_id  = syndromes[(i-2+xSize)%xSize][(j+1)%ySize][t+1]->syndromeID;

          // Edge w. p. 1/2*3/9
          MDist[tc_beta_id][tn_beta_id] = -1*log((double)3/18);
          MDist[tn_beta_id][tc_beta_id] = -1*log((double)3/18);
          if (debug) cout << "New edge " << i << "," << (j+1)%ySize << "," << t+1 << "<->" << i << "," << (j+1)%ySize << "," << t << " wt. " << -1*log((double)3/18) << endl;
          // Edge w. p. 1/2*4/9
          MDist[tc_east_id][tn_beta_id] = -1*log((double)4/18);
          MDist[tn_beta_id][tc_east_id] = -1*log((double)4/18);
          if (debug) cout << "New edge " << i << "," << (j+1)%ySize << "," << t+1 << "<->" << (i+2)%xSize << "," << (j-1+ySize)%ySize << "," << t << " wt. " << -1*log((double)4/18) << endl;
          // Edge w. p. 1/2*5/9
          MDist[tn_beta_id][tn_west_id] = -1*log((double)5/18);
          MDist[tn_west_id][tn_beta_id] = -1*log((double)5/18);
          if (debug) cout << "New edge " << i << "," << (j+1)%ySize << "," << t+1 << "<->" << (i-2+xSize)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)5/18) << endl;
          // Edge w. p. 1/2*6/9
	  MDist[tc_alpha_id][tn_alpha_id] = -1*log((double)6/18);
          MDist[tn_alpha_id][tc_alpha_id] = -1*log((double)6/18);
          if (debug) cout << "New edge " << i << "," << (j-1+ySize)%ySize << "," << t << "<->" << i << "," << (j-1+ySize)%ySize << "," << t+1 << " wt. " << -1*log((double)6/18) << endl;
          // Edge w. p. 1/2*8/9
          MDist[tn_beta_id][tn_alpha_id] = -1*log((double)8/18);
          MDist[tn_alpha_id][tn_beta_id] = -1*log((double)8/18);
          if (debug) cout << "New edge " << i << "," << (j+1)%ySize << "," << t+1 << "<->" << i << "," << (j-1+ySize)%ySize << "," << t+1 << " wt. " << -1*log((double)8/18) << endl;
        }

        // Add edges - circuit model, data
        if ((corLeak == circuit && errType == Z && ((i%2) && !(j%2))) ||
            (corLeak == circuit && errType == X && (!(i%2) && (j%2)))) {
          int tc_1_id  = syndromes[(i-1+xSize)%xSize][j][t]->syndromeID;
          int tc_2_id  = syndromes[(i+1)%xSize][j][t]->syndromeID;
          // Edge w. p. 1/2*3/9
          MDist[tc_1_id][tc_2_id] = -1*log((double)3/18);
          MDist[tc_2_id][tc_1_id] = -1*log((double)3/18);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << j << "," << t << "<->" << (i+1)%xSize << "," << j << "," << t << " wt. " << -1*log((double)3/18) << endl;
          if (t >= TSize - 1)
            break;
          int tn_1_id = syndromes[(i-1+xSize)%xSize][j][t+1]->syndromeID;
          int tn_2_id = syndromes[(i+1)%xSize][j][t+1]->syndromeID;
          // Edge w. p. 1/2*4/9
          MDist[tc_2_id][tn_2_id] = -1*log((double)4/18);
          MDist[tn_2_id][tc_2_id] = -1*log((double)4/18);
          if (debug) cout << "New edge " << (i+1)%xSize << "," << j << "," << t << "<->" << (i+1)%xSize << "," << j << "," << t+1 << " wt. " << -1*log((double)4/18) << endl;
          // Edge w. p. 1/2*5/9
          MDist[tc_1_id][tn_1_id] = -1*log((double)5/18);
          MDist[tn_1_id][tc_1_id] = -1*log((double)5/18);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << j << "," << t << "<->" << (i-1+xSize)%xSize << "," << j << "," << t+1 << " wt. " << -1*log((double)5/18) << endl;
          // Edge w. p. 1/2*8/9
          MDist[tn_1_id][tn_2_id] = -1*log((double)8/18);
          MDist[tn_2_id][tn_1_id] = -1*log((double)8/18);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << j << "," << t+1 << "<->" << (i+1)%xSize << "," << j << "," << t+1 << " wt. " << -1*log((double)8/18) << endl;
        }
      }
    }
  }


  // Update edge disgtances for current leakage history - QUICK MODEL
  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      for (int t = 0; t < TSize; t++) {
        if ((i%2 && !(j%2)) || (!(i%2) && j%2) || !lkgHistory[i][j][t])
          continue;
        if (debug && corLeak == quick) cout << "Updating edges for leakage at i,j,t: " << i << " " << j << " " << t << endl;
        // Add edges - quick model
	if ((corLeak == quick && errType == Z && ((i%2) && (j%2))) || 
	    (corLeak == quick && errType == X && (!(i%2) && !(j%2)))) {
          int tc_3_id = syndromes[(i-1+xSize)%xSize][(j+1)%ySize][t]->syndromeID;
          int tc_4_id = syndromes[(i+1)%xSize][(j+1)%ySize][t]->syndromeID;
          // Edge w. p. 1/2*6/11
          MDist[tc_3_id][tc_4_id] = -1*log((double)6/22);
          MDist[tc_4_id][tc_3_id] = -1*log((double)6/22);
          if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t << " wt. " << -1*log((double)6/22) << endl;
          if (t >= 1) {
            int tp_1_id  = syndromes[(i-1+xSize)%xSize][(j-1+ySize)%ySize][t-1]->syndromeID;
            int tp_2_id  = syndromes[(i+1)%xSize][(j-1+ySize)%ySize][t-1]->syndromeID;
            // Edge w. p. 1/2*1/11
            MDist[tp_1_id][tp_2_id] = -1*log((double)1/22);
            MDist[tp_2_id][tp_1_id] = -1*log((double)1/22);
            if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j-1+ySize)%ySize << "," << t-1 << "<->" << (i+1)%xSize << "," << (j-1+ySize)%ySize << "," << t-1 << " wt. " << -1*log((double)1/22) << endl;
            // Edge w. p. 1/2*2/11
            MDist[tp_1_id][tc_3_id] = -1*log((double)2/22);
            MDist[tc_3_id][tp_1_id] = -1*log((double)2/22);
            if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j-1+ySize)%ySize << "," << t-1 << "<->" << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t << " wt. " << -1*log((double)2/22) << endl;
            // Edge w. p. 1/2*3/11
            MDist[tp_2_id][tc_4_id] = -1*log((double)3/22);
            MDist[tc_4_id][tp_2_id] = -1*log((double)3/22);
            if (debug) cout << "New edge " << (i+1)%xSize << "," << (j-1+ySize)%ySize << "," << t-1 << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t << " wt. " << -1*log((double)3/22) << endl;
	  }
          if (t <= TSize - 2) {
            int tn_3_id = syndromes[(i-1+xSize)%xSize][(j+1)%ySize][t+1]->syndromeID;
            int tn_4_id = syndromes[(i+1)%xSize][(j+1)%ySize][t+1]->syndromeID;
            // Edge w. p. 1/2*7/11
            MDist[tc_4_id][tn_4_id] = -1*log((double)7/22);
            MDist[tn_4_id][tc_4_id] = -1*log((double)7/22);
            if (debug) cout << "New edge " << (i+1)%xSize << "," << (j+1)%ySize << "," << t << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)7/22) << endl;
            // Edge w. p. 1/2*8/11
            MDist[tc_3_id][tn_3_id] = -1*log((double)8/22);
            MDist[tn_3_id][tc_3_id] = -1*log((double)8/22);
            if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t << "<->" << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)8/22) << endl;
            // Edge w. p. 1/2*10/11
            MDist[tn_3_id][tn_4_id] = -1*log((double)10/22);
            MDist[tn_4_id][tn_3_id] = -1*log((double)10/22);
            if (debug) cout << "New edge " << (i-1+xSize)%xSize << "," << (j+1)%ySize << "," << t+1 << "<->" << (i+1)%xSize << "," << (j+1)%ySize << "," << t+1 << " wt. " << -1*log((double)10/22) << endl;
          }
        }

        // Add edges - quick model
	if ((corLeak == quick && errType == Z && (!(i%2) && !(j%2))) ||
            (corLeak == quick && errType == X && ((i%2) && (j%2)))) {
          int tp_alpha_id = 0, tp_beta_id = 0, tp_star_id = 0;
          if (t >= 1) {
            tp_alpha_id       = syndromes[i][(j-2+ySize)%ySize][t-1]->syndromeID;
            tp_beta_id        = syndromes[(i-2+xSize)%xSize][j][t-1]->syndromeID;
            tp_star_id        = syndromes[i][j][t-1]->syndromeID;
          }
          int tc_star_id       = syndromes[i][j][t]->syndromeID;
          int tc_gamma_id      = syndromes[(i+2)%xSize][j][t]->syndromeID;
          int tc_delta_id      = syndromes[i][(j+2)%ySize][t]->syndromeID;
          int tn_delta_id = 0, tn_star_id = 0, tn_deltawest_id = 0;
          if (t <= TSize - 2) {
            tn_delta_id      = syndromes[i][(j+2)%ySize][t+1]->syndromeID;
            tn_star_id       = syndromes[i][j][t+1]->syndromeID;
            tn_deltawest_id  = syndromes[(i-2+xSize)%xSize][(j+2)%ySize][t+1]->syndromeID;
          }
          if (t >= 1 && t <= TSize - 2) {
            // Edge w. p. 1/2*1/11
            MDist[tp_alpha_id][tn_delta_id] = -1*log((double)1/22);
            MDist[tn_delta_id][tp_alpha_id] = -1*log((double)1/22);
            if (debug) cout << "New edge " << i << "," << (j-2+ySize)%ySize << "," << t-1 << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)1/22) << endl;
            // Edge w. p. 1/2*2/11
            MDist[tp_beta_id][tn_delta_id] = -1*log((double)2/22);
            MDist[tn_delta_id][tp_beta_id] = -1*log((double)2/22);
            if (debug) cout << "New edge " << (i-2+xSize)%xSize << "," << j << "," << t-1 << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)2/22) << endl;
          }
          if (t >= 1 && t <= TSize - 1) {
            // Edge w. p. 1/2*5/11
            MDist[tp_star_id][tc_star_id] = -1*log((double)5/22);
            MDist[tc_star_id][tp_star_id] = -1*log((double)5/22);
            if (debug) cout << "New edge " << i << "," << j << "," << t-1 << "<->" << i << "," << j << "," << t+1 << " wt. " << -1*log((double)5/22) << endl;
          }
          if (t <= TSize - 2) {
            // Edge w. p. 1/2*6/11
            MDist[tc_delta_id][tn_delta_id] = -1*log((double)6/22);
            MDist[tn_delta_id][tc_delta_id] = -1*log((double)6/22);
            if (debug) cout << "New edge " << i << "," << (j+2)%ySize << "," << t << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)6/22) << endl;
            // Edge w. p. 1/2*(7+3)/11 - 3/22*7/22
            MDist[tc_gamma_id][tn_delta_id] = -1*log((double)199/484);
            MDist[tn_delta_id][tc_gamma_id] = -1*log((double)199/484);
            if (debug) cout << "New edge " << (i+2)%xSize << "," << j << "," << t << "<->" << i << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)199/484) << endl;
            // Edge w. p. 1/2*8/11
            MDist[tn_delta_id][tn_deltawest_id] = -1*log((double)8/22);
            MDist[tn_deltawest_id][tn_delta_id] = -1*log((double)8/22);
            if (debug) cout << "New edge " << i << "," << (j+2)%ySize << "," << t+1 << "<->" << (i-2+xSize)%xSize << "," << (j+2)%ySize << "," << t+1 << " wt. " << -1*log((double)8/22) << endl;
            // Edge w. p. 1/2*10/11
            MDist[tn_delta_id][tc_star_id] = -1*log((double)10/22);
            MDist[tc_star_id][tn_delta_id] = -1*log((double)10/22);
            if (debug) cout << "New edgee " << i << "," << (j+2)%ySize << "," << t+1 << "<->" << i << "," << j << "," << t+1 << " wt. " << -1*log((double)10/22) << endl;

          }
        }
      }
    }
  }

  // Perform Floyd Warshall
  for (int k = 0; k < numSynd; k++) {
    for (int i = 0; i < numSynd; i++) {
      for (int j = 0; j < numSynd; j++) {
        if (MDist[i][k] + MDist[k][j] < MDist[i][j])
	  MDist[i][j] = MDist[i][k] + MDist[k][j];
      }      
    }
  }

  // For each vertex u, connect to every other vertex v by an edge. No edge prunning performed.
  for (int u = 0; u < (int)vertices.size(); u++) {
    for (int v = 0; v < (int)vertices.size(); v++) {
      if (u == v)
        continue;
      // Add the edge with its distance
      int uid = vertices[u]->id;
      int vid = vertices[v]->id;
      int usid = vertices[u]->syndromeID;
      int vsid = vertices[v]->syndromeID;
      double dist = MDist[usid][vsid];
      MyEdge* newEdge = new MyEdge (uid,vid,dist);
      edges.push_back(newEdge);
      if (debug) {
        int ux  = vertices[u]->x;
        int uy  = vertices[u]->y;
        int ut  = vertices[u]->t;
        int vx  = vertices[v]->x;
        int vy  = vertices[v]->y;
        int vt  = vertices[v]->t;
        cout << "Adding new edge: uid=" << uid << "(" << ux << "," << uy << "," << ut << ") and vid=" << vid << "(" << vx << "," << vy << "," << vt << "): dist=" << dist <<  endl;
      }

    }
  }
}


// Calculate distance of two syndromes in the grid with edges of type a, b, c, d, and e.
double Cube::getDistance (int x1Loc, int y1Loc, int t1Loc, int x2Loc, int y2Loc, int t2Loc) {

  double dist = 0;
  int delta = 0;

  // Initialize x1, y1, t1,  x2, y2, t2 so that t1 > t2
  int x1 = x1Loc, y1 = y1Loc, t1 = t1Loc, x2 = x2Loc, y2 = y2Loc, t2 = t2Loc;
  if (t2Loc > t1Loc) {
    x1 = x2Loc; y1 = y2Loc; t1 = t2Loc; x2 = x1Loc; y2 = y1Loc; t2 = t1Loc;
  }
  if (debug) cout << "getDistance(x1=" << x1 << " y1=" << y1 << " t1=" << t1 << " x2=" << x2 << " y2=" << y2 << " t2=" << t2 << ")";
 
  // While NE move desired and t1-t2>0 use f edge
  delta = min_3int (t1-t2, min(abs(x1-x2), xSize-abs(x1-x2)) / 2, min(abs(y1-y2), ySize-abs(y1-y2)) / 2);
  if ( delta &&
       ((x1 < x2 && x2 - x1 <= xSize - x2 + x1) ||  // east direction
        (x2 < x1 && x1 - x2 >= xSize - x1 + x2)) &&
       ((y1 < y2 && y2 - y1 >= ySize - y2 + y1) ||  // north direction
        (y2 < y1 && y1 - y2 <= ySize - y1 + y2))) {
    x1 = (x1 + 2*delta) % (xSize);
    y1 = (y1 - 2*delta + ySize) % (ySize);
    t1 -= delta; 
    dist += delta * wf;
    if (debug) cout << " uses edges " << delta << "*f";
  }

  // While N move desired and t1-t2>0 use c edge
  delta = min (t1-t2, min(abs(y1-y2), ySize-abs(y1-y2)) / 2 );
  if ( delta &&
       ((y1 < y2 && y2 - y1 >= ySize - y2 + y1) ||  // north direction
        (y2 < y1 && y1 - y2 <= ySize - y1 + y2))) {
    y1 = (y1 - 2*delta + ySize) % (ySize);
    t1 -= delta; 
    dist += delta * wc;
    if (debug) cout << " " << delta << "*c";
  }

  // While W move desired and t1-t2>0 use e edge
  delta = min (t1-t2, min(abs(x1-x2), xSize-abs(x1-x2)) / 2 );
  if ( delta &&
       ((x1 < x2 && x2 - x1 >= xSize - x2 + x1) ||  // west direction              
	(x2 < x1 && x1 - x2 <= xSize - x1 + x2))) {
    x1 = (x1 - 2*delta + xSize) % (xSize);
    t1 -= delta;
    dist += delta * we;
    if (debug) cout << " " << delta << "*e";
  }

  // While t1-t2>0 use a edge
  delta = t1 - t2;
  if ( delta ) {
    dist += delta * wa;
    if (debug) cout << " " << delta << "*a";
  }

  // While N or S move desired use b edge
  delta = min(abs(y1-y2), ySize-abs(y1-y2)) / 2;
  if ( delta ) {
    dist += delta * wb;
    if (debug) cout << " " << delta << "*b";
  }

  // While W or E move desired use d edge
  delta = min(abs(x1-x2), xSize-abs(x1-x2)) / 2;
  if ( delta ) {
    dist += delta * wd;
    if (debug) cout << " " << delta << "*d";
  }

  if (debug) cout << ": total distance " << dist << endl;
  return dist;

}


// Return the smallest element in a vector of doubles.
double Cube::min_double (vector<double> & vd) {

  assert (vd.size() > 0);

  double result = vd[0];
  for (int i = 0; i < (int)vd.size(); i++) {
    result = min (result, vd[i]);
  }

  return result;

}


// Return the smallest of three integers.
int Cube::min_3int (int x, int y, int z) {

  int m = min (x, y);
  return min (m, z);

}


// Constructor generates a unit cell. In case of the toric code this is a single qubit.
Cell::Cell (void) {
  ErrCurrent = new Operator (1);
}


// Destructor for a unit cell.
Cell::~Cell (void) {
  delete ErrCurrent;
}


/*! Constructor generates the lattice.
 * \param newXsize minimum distance of the toric code
 * \param newT number of syndrome measurements
 * \param pError error probability of idle timestep
 * \param p2Error error probability of two qubit gate
 * \param qError error probabilty of measurement
 * \param pLeakUp probability of erasure
 * \param pLeakDown probability of "unerasure"
 * \param circuitErr if 0 use phenomenological model, otherwise use circuit model
 */
Lattice::Lattice(int newXSize, int newT, double pError, double p2Error, double qError, double pLeakUp, double pLeakDown, double addlZ, bool circuitErr, corTech correctLeak, Counter *pCtr) {

  xSize = 2*newXSize;
  ySize = 2*newXSize;
  T = newT;
  p = pError;
  p2 = p2Error;
  q = qError;
  pu = pLeakUp;
  pd = pLeakDown;
  addZ = addlZ;
  circErr = circuitErr;
  corLeak = correctLeak;
  ctr = pCtr;
  zero = 0.0;
 
  cells = new Cell**[xSize];
  for (int i = 0; i < xSize; i++) {
    cells[i] = new Cell*[ySize];
    for (int j = 0; j < ySize; j++) {
      // Space is also reserved for ancillas. If no ancillas desired uncomment.
      //if ((i%2 && j%2) || (!(i%2) && !(j%2)))
      //  continue;
      cells[i][j] = new Cell ();
    }
  }

  // Add leakage -- corresponding to the equilibrium
  double eqRateData = 0;
  double eqRateAnc  = 0;
  if (pu > 0 && corLeak == no) {
    /* It is possible to calculate the equilibrium iteratively
    double pl = 0, plPrev = -1;
    while (pl - plPrev > 0.001*pu) {
      plPrev = pl;
      pl = pl*(1-pd);
      pl = (pl+(1-pl)*pu)*(1-pd);
      pl = (pl+(1-pl)*pu)*(1-pd);
      pl = (pl+(1-pl)*pu)*(1-pd);
      pl = (pl+(1-pl)*pu)*(1-pd);
      pl = pl*(1-pd);
    }*/
    // Here we calculate the equilibrium up to the highest order terms
    eqRateData = 4*pu/(4*pu + 6*pd);
    eqRateAnc  = 0;
  }
  if (pu > 0 && corLeak == quick) {
    double pl = 0;
    for (int i = 0; i < 5; i++) {
      pl = (pl+(1-pl)*pu)*(1-pd);
    } 
    eqRateData = pl;
    eqRateAnc  = 0;
  }
  if (pu > 0 && (corLeak == gate || corLeak == circuit)) {
    double pl = 0;
    for (int i = 0; i < 3; i++) {
      pl = (pl+(1-pl)*pu)*(1-pd);
    }
    eqRateData = pl;
    eqRateAnc  = 0;
  }
  if (debug) cout << "eqRateData=" << eqRateData << " eqRateAnc=" << eqRateAnc << endl;
  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      if (!ctr->usePath && !ctr->countPath) {
	double r = myRandom();
        // If ancilla...
        if (((i%2 && j%2) || (!(i%2) && !(j%2))) && r < eqRateAnc)
          cells[i][j]->ErrCurrent->ops[0] = L;
        // If data...
        if (!((i%2 && j%2) || (!(i%2) && !(j%2))) && r < eqRateData)
          cells[i][j]->ErrCurrent->ops[0] = L;
      }
    }
  }

}


// Destructor of the lattice.
Lattice::~Lattice(void) {

  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      //if ((i%2 && j%2) || (!(i%2) && !(j%2)))
      //  continue;
      delete cells[i][j];
    }
    delete[] cells[i];
  }
  delete[] cells;

}


// Wrapper for the simulation
void Lattice::simulate (){

  // Simulates T rounds of error propagation and one round of perfect measurement, finished by the 
  // error correction step. Calls the mathcing algorithm to find pairs of squares with nontrivial 
  // syndromes to be corrected.
  if (debug) cout << "SIMULATING T ROUNDS OF ERROR PROPAGATION:" << endl;

  // The cube maintains syndrome history throughout the T + 1 time steps
  Cube cube (xSize, ySize, T + 1);
  cube.corLeak = corLeak;

  // Simulate the T step of error propagations and measurements
  for (int t = 0; t < T; t++) {
    
    // Introduce new error - phenomenological model only
    if (!circErr) {
      for (int i = 0; i < xSize; i++) {
        for (int j = 0; j < ySize; j++) {
          if ((i%2 && j%2) || (!(i%2) && !(j%2)))
            continue;
          cells[i][j]->ErrCurrent->addSingleDepol (p, addZ, 0, ctr);
        }
      }
      // Print the current error state
      if (debug) cout << "ERRORS AT TIME t=" << t << ":" << endl;
      if (debug) printState();
    }

    // Obtain syndromes - phenomenological model only
    if (!circErr) {
      for (int i = 0; i < xSize; i++) {
        for (int j = 0; j < ySize; j++) {
          if ( i%2 && j%2 && phenomSyndromeZZZZ(i, j, q) )
            cube.syndromes[i][j][t]->value = true;
          if ( !(i%2) && !(j%2) && phenomSyndromeXXXX(i, j, q) )
            cube.syndromes[i][j][t]->value = true;
        }
      }
    }

    // Obtain syndromes - circuit model
    if (circErr) {
      // Step 1 (initialization)
      for (int i = 0; i < xSize; i++) {
        for (int j = 0; j < ySize; j++) {
          if ( i%2 && j%2 )
            circSyndromeZZZZ_1(i, j, p, pu, pd);
          if ( !(i%2) && !(j%2) )
            circSyndromeXXXX_1(i, j, p, pu, pd);
        }
      }
      if (debug) cout << "CURRENT STATE t=" << t << " step=1:" << endl;
      if (debug) printState();

      // Step 2 (1st CNOT)
      for (int i = 0; i < xSize; i++) {
        for (int j = 0; j < ySize; j++) {
          if ( i%2 && j%2 )
            circSyndromeZZZZ_2(i, j, p, pu, pd);
          if ( !(i%2) && !(j%2) )
            circSyndromeXXXX_2(i, j, p, pu, pd);
        }
      }
      if (debug) cout << "CURRENT STATE t=" << t << " step=2:" << endl;
      if (debug) printState();

      // Step 3 (2nd CNOT)
      for (int i = 0; i < xSize; i++) {
        for (int j = 0; j < ySize; j++) {
          if ( i%2 && j%2 )
            circSyndromeZZZZ_3(i, j, p, pu, pd);
          if ( !(i%2) && !(j%2) )
            circSyndromeXXXX_3(i, j, p, pu, pd);
        }
      }
      if (debug) cout << "CURRENT STATE t=" << t << " step=3:" << endl;
      if (debug) printState();

      // Step 4 (3rd CNOT)
      for (int i = 0; i < xSize; i++) {
        for (int j = 0; j < ySize; j++) {
          if ( i%2 && j%2 )
            circSyndromeZZZZ_4(i, j, p, pu, pd);
          if ( !(i%2) && !(j%2) )
            circSyndromeXXXX_4(i, j, p, pu, pd);
        }
      }
      if (debug) cout << "CURRENT STATE t=" << t << " step=4:" << endl;
      if (debug) printState();

      // Step 5 (4th CNOT)
      for (int i = 0; i < xSize; i++) {
        for (int j = 0; j < ySize; j++) {
          if ( i%2 && j%2 )
            circSyndromeZZZZ_5(i, j, p, pu, pd);
          if ( !(i%2) && !(j%2) )
            circSyndromeXXXX_5(i, j, p, pu, pd);
        }
      }
      if (debug) cout << "CURRENT STATE t=" << t << " step=5:" << endl;
      if (debug) printState();

      // Step 6 (measurement)
      int XXXXCnt = 0, ZZZZCnt = 0;
      bool lkgFlagA = false, lkgFlagD = false;
      for (int i = 0; i < xSize; i++) {
        for (int j = 0; j < ySize; j++) {
          if (i%2 && j%2) {
            if (circSyndromeZZZZ_6(i, j, p, q, pu, pd, &lkgFlagA, &lkgFlagD)) {
              cube.syndromes[i][j][t]->value = true;
              ZZZZCnt++;
            }
            cube.lkgHistory[i][j][t] = lkgFlagA;
            cube.lkgHistory[i][(j+1)%ySize][t] = lkgFlagD;
          }
          if (!(i%2) && !(j%2)) {
            if (circSyndromeXXXX_6(i, j, p, q, pu, pd, &lkgFlagA, &lkgFlagD)) {
              cube.syndromes[i][j][t]->value = true;
              XXXXCnt++;
            }
            cube.lkgHistory[i][j][t] = lkgFlagA;
            cube.lkgHistory[i][(j+1)%ySize][t] = lkgFlagD;
          }
        }
      }
      if (debug) cout << "CURRENT STATE t=" << t << " step=6:" << endl;
      if (debug) printState();
      if (debug) cout << "Number of measured syndromes: XXXXCnt=" << XXXXCnt << " ZZZZCnt=" << ZZZZCnt << endl;
    }

    // Summarize current errors and syndromes
    log << "state after round at time " << t << endl;
    for (int j = 0; j < ySize; j++) {
      for (int i = 0; i < xSize; i++) {
        if ((i%2&&j%2)||(!(i%2)&&!(j%2)))  
          log << cube.syndromes[i][j][t]->value;
        else
          switch(cells[i][j]->ErrCurrent->ops[0]) {
            case I: log << "."; break;
            case X: log << "X"; break;
            case Y: log << "Y"; break;
            case Z: log << "Z"; break;
            case L: log << "L"; break;
            default: assert(false);
          }
      }
      log << endl;
    }

  }

  // If this is the last round of the simulation, remove leakage,
  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      if ((i%2 && j%2) || (!(i%2) && !(j%2)))
        continue;
      cells[i][j]->ErrCurrent->removeLeakage();
    }
  }
  // and add one last round of perfect measurement
  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      if ( i%2 && j%2 && phenomSyndromeZZZZ(i, j, 0, false) )
        cube.syndromes[i][j][T]->value = true;
      if ( !(i%2) && !(j%2) && phenomSyndromeXXXX(i, j, 0, false) )
        cube.syndromes[i][j][T]->value = true;
    }
  }

  // Summarize current errors and syndromes
  log << "state after ideal round at time " << T << endl;
  for (int j = 0; j < ySize; j++) {
    for (int i = 0; i < xSize; i++) {
      if ((i%2&&j%2)||(!(i%2)&&!(j%2)))  
        log << cube.syndromes[i][j][T]->value;
      else
        switch(cells[i][j]->ErrCurrent->ops[0]) {
          case I: log << "."; break;
          case X: log << "X"; break;
          case Y: log << "Y"; break;
          case Z: log << "Z"; break;
          case L: log << "L"; break;
          default: assert(false);
        }
    }
    log << endl;
  }

  // Calculate relative edge weights of time edges, space edges, and diagonals
  if (ctr->usePath)
    cube.accumulateRelativeWeights (circErr, ctr);

  if (!ctr->usePath && !ctr->countPath)
    cube.finalizeRelativeWeights (p, p2, q, circErr, ctr);

  // Create prunned vertex and edge list
  vector<Vertex*> verticesXXXX;
  vector<Vertex*> verticesZZZZ;
  vector<MyEdge*> edgesXXXX;
  vector<MyEdge*> edgesZZZZ;
  cube.getVertices(verticesXXXX, verticesZZZZ);
  cube.initNearest();
  if (debug) cout << "OBTAINING XXXX EDGES:" << endl;
  if (lkgDetect && !ctr->usePath && !ctr->countPath)
    cube.getEdgesLkgDetBoost(verticesXXXX, edgesXXXX, Z);
  else
    cube.getEdges(verticesXXXX, edgesXXXX);
  if (debug) cout << "OBTAINING ZZZZ EDGES:" << endl;
  if (lkgDetect && !ctr->usePath && !ctr->countPath)
    cube.getEdgesLkgDetBoost(verticesZZZZ, edgesZZZZ, X);
  else
    cube.getEdges(verticesZZZZ, edgesZZZZ);

  // Call the matching subroutine that also corrects the detected errors
  if (debug) cout << "MATCHING ON ZZZZ SYNDROMES:" << endl;
  callMatching (X, verticesZZZZ, edgesZZZZ);
  if (debug) cout << "MATCHING ON XXXX SYNDROMES:" << endl;
  callMatching (Z, verticesXXXX, edgesXXXX);

  // Summarize final errors after correction
  log << "state after corrections" << endl;
  for (int j = 0; j < ySize; j++) {
    for (int i = 0; i < xSize; i++) {
      if ((i%2&&j%2)||(!(i%2)&&!(j%2)))  
        log << " ";
      else
        switch(cells[i][j]->ErrCurrent->ops[0]) {
          case I: log << "."; break;
          case X: log << "X"; break;
          case Y: log << "Y"; break;
          case Z: log << "Z"; break;
          case L: log << "L"; break;
          default: assert(false);
        }
    }
    log << endl;
  }

  // Delete edges and vertices
  for (int i = 0; i < (int) verticesZZZZ.size(); i++) {
    delete verticesZZZZ[i];
  }
  for (int i = 0; i < (int) verticesXXXX.size(); i++) {
    delete verticesXXXX[i];
  }
  for (int i = 0; i < (int) edgesZZZZ.size(); i++) {
    delete edgesZZZZ[i];
  }
  for (int i = 0; i < (int) edgesXXXX.size(); i++) {
    delete edgesXXXX[i];
  }

}


// Calculates site syndrome at (around) the specified coordinates. Phenomenological model.
bool Lattice::phenomSyndromeXXXX (int xLoc, int yLoc, double qProb, bool noisy){

  assert (xLoc < xSize && yLoc < ySize && !(xLoc % 2) && !(yLoc % 2));

  bool result = false;

  Operator *S = new Operator(4);
  S->ops[0]  = X;
  S->ops[1]  = X;
  S->ops[2]  = X;
  S->ops[3]  = X;

  Operator *site = new Operator(4);
  site->ops[0]  = cells[xLoc][(yLoc-1+ySize)%ySize]->ErrCurrent->ops[0];
  site->ops[1]  = cells[(xLoc+1)%xSize][yLoc]->ErrCurrent->ops[0];
  site->ops[2]  = cells[xLoc][(yLoc+1)%ySize]->ErrCurrent->ops[0];
  site->ops[3]  = cells[(xLoc-1+xSize)%xSize][yLoc]->ErrCurrent->ops[0];

  result = !S->commute(*site);

  if(noisy) {
    if(ctr->usePath) {
      if (ctr->modeMap[ctr->locCtr] == 1) {
        result = !result;
        if (debug) cout << "Syndrome measurement error occurs: XXXX at x=" << xLoc << " y=" << yLoc << endl;
      }
    } else if (!ctr->usePath && !ctr->countPath) {
      // The syndrome is incorrect with probability q.
      double r = myRandom();
      if (r < qProb) {
        result = !result;
        if (debug) cout << "Syndrome measurement error occurs: XXXX at x=" << xLoc << " y=" << yLoc << endl;
      }
    }
    ctr->locModes.push_back(2);
    ctr->locCauses.push_back(PhenSynd);
    ctr->locProbs.push_back(qProb);
    ctr->locCtr++;
  }

  delete S;
  delete site;
  return result;
}


// Calculates plaquette syndrome at (around) the specified coordinates. Phenomenological model.
bool Lattice::phenomSyndromeZZZZ (int xLoc, int yLoc, double qProb, bool noisy){

  assert (xLoc < xSize && yLoc < ySize && (xLoc % 2) && (yLoc % 2));

  bool result = false;

  Operator *P = new Operator(4);
  P->ops[0]  = Z;
  P->ops[1]  = Z;
  P->ops[2]  = Z;
  P->ops[3]  = Z;

  Operator *plaquette = new Operator(4);
  plaquette->ops[0]  = cells[xLoc][yLoc-1]->ErrCurrent->ops[0];
  plaquette->ops[1]  = cells[(xLoc+1)%xSize][yLoc]->ErrCurrent->ops[0];
  plaquette->ops[2]  = cells[xLoc][(yLoc+1)%ySize]->ErrCurrent->ops[0];
  plaquette->ops[3]  = cells[xLoc-1][yLoc]->ErrCurrent->ops[0];

  result = !P->commute(*plaquette);

  if(noisy) {
    if(ctr->usePath) {
      if (ctr->modeMap[ctr->locCtr] == 1) {
        result = !result;
        if (debug) cout << "Syndrome measurement error occurs: ZZZZ at x=" << xLoc << " y=" << yLoc << endl;
      }
    } else if (!ctr->usePath && !ctr->countPath) {
      // The syndrome is incorrect with probability q.
      double r = myRandom();
      if (r < qProb) {
        result = !result;
        if (debug) cout << "Syndrome measurement error occurs: ZZZZ at x=" << xLoc << " y=" << yLoc << endl;
      }
    }
    ctr->locModes.push_back(2);
    ctr->locCauses.push_back(PhenSynd);
    ctr->locProbs.push_back(qProb);
    ctr->locCtr++;
  }

  delete P;
  delete plaquette;
  return result;
}


// Leakage reduction circuit that initializes a new unleaked qubit and swaps it
// with the current qubit.
void Lattice::lkgRedCircuit (Cell* data, double pProb, double puProb, double pdProb, bool *lkgFlagD) {
  Cell* dataTmp = new Cell ();
  // Swap
  data->ErrCurrent->SWAP (pProb, puProb, pdProb, *dataTmp->ErrCurrent, 0, 0, ctr);
  // Detect leakage
  if (data->ErrCurrent->ops[0] == L && lkgFlagD != NULL)
    *lkgFlagD = true;
  else if (data->ErrCurrent->ops[0] != L && lkgFlagD != NULL)
    *lkgFlagD = false;
  // Move back
  data->ErrCurrent->ops[0] = dataTmp->ErrCurrent->ops[0];
  delete dataTmp;
}


// Calculates site syndrome at (around) the specified coordinates. Simulates leakage.
void Lattice::circSyndromeXXXX_1 (int xLoc, int yLoc, double pProb, double puProb, double pdProb){

  assert (xLoc < xSize && yLoc < ySize && !(xLoc % 2) && !(yLoc % 2));

  // Locate the qubit in the up direction and the ancilla
  Cell *data1   = cells[xLoc][(yLoc-1+ySize)%ySize];
  Cell *ancilla = cells[xLoc][yLoc];

  // Re-initialize the ancilla, Z error probability is p
  ancilla->ErrCurrent->initQubit (pProb, zero, 0, Z, ctr);

  // Simulate leakage of data
  data1->ErrCurrent->identityErrLkg (pProb, zero, zero, addZ, 0, ctr);

}

void Lattice::circSyndromeXXXX_2 (int xLoc, int yLoc, double pProb, double puProb, double pdProb){

  // Locate the first data qubit (up) and the ancilla
  Cell *data1   = cells[xLoc][(yLoc-1+ySize)%ySize];
  Cell *ancilla = cells[xLoc][yLoc];

  // Simulate error propagation from the first CNOT, inluding leakage of data and ancilla
  ancilla->ErrCurrent->CNOTErrLkg (pProb, zero, zero, *data1->ErrCurrent, 0, 0, ctr);
    
  // Simulate Magnetic Field Flucuations
  data1->ErrCurrent->identityErrLkg (pProb, zero, zero, addZ, 0, ctr);

  // Gate leakage reduction
  if (corLeak == gate) {
    lkgRedCircuit (data1, pProb, puProb, pdProb, NULL);
    lkgRedCircuit (ancilla, pProb, puProb, pdProb, NULL);
  }
}

void Lattice::circSyndromeXXXX_3 (int xLoc, int yLoc, double pProb, double puProb, double pdProb){

  // Locate the decond data qubit (left) and the ancilla
  Cell *data2   = cells[(xLoc-1+xSize)%xSize][yLoc];
  Cell *ancilla = cells[xLoc][yLoc];

  // Simulate error propagation from the second CNOT, inluding leakage of data and ancilla
  ancilla->ErrCurrent->CNOTErrLkg (pProb, puProb, pdProb, *data2->ErrCurrent, 0, 0, ctr);
    
  // Simulate Magnetic Field Flucuations
  data2->ErrCurrent->identityErrZ (pProb, addZ, 0, ctr);
    


  // Gate leakage reduction
  if (corLeak == gate) {
    lkgRedCircuit (data2, pProb, puProb, pdProb, NULL);
    lkgRedCircuit (ancilla, pProb, puProb, pdProb, NULL);
  }
}

void Lattice::circSyndromeXXXX_4 (int xLoc, int yLoc, double pProb, double puProb, double pdProb){

  // Locate the third data qubit (right) and the ancilla
  Cell *data3   = cells[(xLoc+1)%xSize][yLoc];
  Cell *ancilla = cells[xLoc][yLoc];

  // Simulate error propagation from the third CNOT, inluding leakage of data and ancilla
  ancilla->ErrCurrent->CNOTErrLkg (pProb, puProb, pdProb, *data3->ErrCurrent, 0, 0, ctr);
    
  // Simulate Magnetic Field Flucuations
  data3->ErrCurrent->identityErrZ (pProb, addZ, 0, ctr);
    


  // Gate leakage reduction
  if (corLeak == gate) {
    lkgRedCircuit (data3, pProb, puProb, pdProb, NULL);
    lkgRedCircuit (ancilla, pProb, puProb, pdProb, NULL);
  }
}

void Lattice::circSyndromeXXXX_5 (int xLoc, int yLoc, double pProb, double puProb, double pdProb){

  // Locate the fourth data qubit (down) and the ancilla
  Cell *data4   = cells[xLoc][(yLoc+1)%ySize];
  Cell *ancilla = cells[xLoc][yLoc];

  if (corLeak == quick) {
    data4->ErrCurrent->CNOTErrLkg (pProb, puProb, pdProb, *ancilla->ErrCurrent, 0, 0, ctr);
  }
  // Simulate error propagation from the fourth CNOT, inluding leakage of data and ancilla
  ancilla->ErrCurrent->CNOTErrLkg (pProb, puProb, pdProb, *data4->ErrCurrent, 0, 0, ctr);
    
    
  // Simulate Magnetic Field Flucuations
  data4->ErrCurrent->identityErrZ (pProb, addZ, 0, ctr);
    


  // Gate leakage reduction only needed on data after the last CNOT
  if (corLeak == gate) {
    lkgRedCircuit (data4, pProb, puProb, pdProb, NULL);
  }

}

bool Lattice::circSyndromeXXXX_6 (int xLoc, int yLoc, double pProb, double qProb, double puProb, double pdProb, bool *lkgFlagA, bool *lkgFlagD){

  // Locate the data qubit in the down direction and the ancilla
  Cell *data4, *ancilla;
  if (corLeak == quick) {
    data4   = cells[xLoc][yLoc];
    ancilla = cells[xLoc][(yLoc+1)%ySize];
  } else {
    data4   = cells[xLoc][(yLoc+1)%ySize];
    ancilla = cells[xLoc][yLoc];
  }
  bool result = false;

  // Obtain syndrome by reading the acnilla
  if (ancilla->ErrCurrent->ops[0] == Y || ancilla->ErrCurrent->ops[0] == Z)
    result = true;

  // If the ancilla is leaked, the syndrome is also true
  if (!lkgDetect && ancilla->ErrCurrent->ops[0] == L) {
    result = true;
  }

  // If leakage is detected
  *lkgFlagA = (ancilla->ErrCurrent->ops[0] == L);
  if (lkgDetect && ancilla->ErrCurrent->ops[0] == L) {
    result = false;
  }

  // Meanwhile data leaks again...
  data4->ErrCurrent->identityErrLkg (pProb, zero, zero, addZ, 0, ctr);

  // The data needs to be moved back before next round of simulation
  if (corLeak == quick)
    ancilla->ErrCurrent->ops[0] = data4->ErrCurrent->ops[0];

  // Here we do leakage reduction for the circuit
  if (corLeak == circuit) {
    lkgRedCircuit (data4, pProb, puProb, pdProb, lkgFlagD);
  }

  // The syndrome is incorrect with probability p
  if (debug) cout << "Calculating XXXX syndrome at xLoc=" << xLoc << " yLoc=" << yLoc << ":"<< endl;
  if (debug) {cout << "  ancilla:"; ancilla->ErrCurrent->printState();}
  if (!lkgDetect || ancilla->ErrCurrent->ops[0] != L) {
    if(ctr->usePath) {
      if (ctr->modeMap[ctr->locCtr] == 1) {
        result = !result;
        if (debug) cout << "syndrome error due to measurement occured" << endl;
      }
    } else if (!ctr->usePath && !ctr->countPath) {
      // The syndrome is incorrect with probability q.                                                       
      double r = myRandom();
      if (r < qProb) {
        result = !result;
        if (debug) cout << "  syndrome error due to measurement occured" << endl;
      }
    }
    ctr->locModes.push_back(2);
    ctr->locCauses.push_back(MeasErr);
    ctr->locProbs.push_back(qProb);
    ctr->locCtr++;
  }
  if (debug) cout << "  syndrome=" << result << endl;

  return result;
}


// Calculates site syndrome at (around) the specified coordinates. Simulates leakage.
void Lattice::circSyndromeZZZZ_1 (int xLoc, int yLoc, double pProb, double puProb, double pdProb){

  assert (xLoc < xSize && yLoc < ySize && (xLoc % 2) && (yLoc % 2));

  // Locate the data qubit in the up direction and the ancilla
  Cell *data1  = cells[xLoc][yLoc-1];
  Cell *ancilla = cells[xLoc][yLoc];

  // Re-initialize the ancilla, X error probability is p
  ancilla->ErrCurrent->initQubit (pProb, puProb, 0, Z, ctr);

  // Simulate leakage of data
  data1->ErrCurrent->identityErrLkg (pProb, zero, zero, addZ, 0, ctr);

}

void Lattice::circSyndromeZZZZ_2 (int xLoc, int yLoc, double pProb, double puProb, double pdProb){

  // Locate the first data qubit (up) and the ancilla
  Cell *data1  = cells[xLoc][yLoc-1];
  Cell *ancilla = cells[xLoc][yLoc];

  // Simulate error propagation from the first CNOT, inluding leakage of data and ancilla
  data1->ErrCurrent->CNOTErrLkg (pProb, puProb, pdProb, *ancilla->ErrCurrent, 0, 0, ctr);
    
    
  // Simulate Magnetic Field Flucuations
  data1->ErrCurrent->identityErrZ (pProb, addZ, 0, ctr);
    


  // Gate leakage reduction
  if (corLeak == gate) {
    lkgRedCircuit (data1, pProb, puProb, pdProb, NULL);
    lkgRedCircuit (ancilla, pProb, puProb, pdProb, NULL);
  }
}

void Lattice::circSyndromeZZZZ_3 (int xLoc, int yLoc, double pProb, double puProb, double pdProb){

  // Locate the second data qubit (left) and the ancilla
  Cell *data2  = cells[xLoc-1][yLoc];
  Cell *ancilla = cells[xLoc][yLoc];

  // Simulate error propagation from the second CNOT, inluding leakage of data and ancilla
  data2->ErrCurrent->CNOTErrLkg (pProb, puProb, pdProb, *ancilla->ErrCurrent, 0, 0, ctr);
    
    
  // Simulate Magnetic Field Flucuations
  data2->ErrCurrent->identityErrZ (pProb, addZ, 0, ctr);
    


  // Gate leakage reduction
  if (corLeak == gate) {
    lkgRedCircuit (data2, pProb, puProb, pdProb, NULL);
    lkgRedCircuit (ancilla, pProb, puProb, pdProb, NULL);
  }
}

void Lattice::circSyndromeZZZZ_4 (int xLoc, int yLoc, double pProb, double puProb, double pdProb){

  // Locate the third data qubit (right) and the ancilla
  Cell *data3  = cells[(xLoc+1)%xSize][yLoc];
  Cell *ancilla = cells[xLoc][yLoc];

  // Simulate error propagation from the third CNOT, inluding leakage of data and ancilla
  data3->ErrCurrent->CNOTErrLkg (pProb, puProb, pdProb, *ancilla->ErrCurrent, 0, 0, ctr);
    
  // Simulate Magnetic Field Flucuations
  data3->ErrCurrent->identityErrZ (pProb, addZ, 0, ctr);
    


  // Gate leakage reduction
  if (corLeak == gate) {
    lkgRedCircuit (data3, pProb, puProb, pdProb, NULL);
    lkgRedCircuit (ancilla, pProb, puProb, pdProb, NULL);
  }
}

void Lattice::circSyndromeZZZZ_5 (int xLoc, int yLoc, double pProb, double puProb, double pdProb){

  // Locate the fourth data qubit (down) and the ancilla
  Cell *data4  = cells[xLoc][(yLoc+1)%ySize];
  Cell *ancilla = cells[xLoc][yLoc];
  
  if (corLeak == quick) {
    ancilla->ErrCurrent->CNOTErrLkg (pProb, puProb, pdProb, *data4->ErrCurrent, 0, 0, ctr);
  }


  // Simulate error propagation from the fourth CNOT, inluding leakage of data and ancilla
  data4->ErrCurrent->CNOTErrLkg (pProb, puProb, pdProb, *ancilla->ErrCurrent, 0, 0, ctr);
    
  // Simulate Magnetic Field Flucuations
  data4->ErrCurrent->identityErrZ (pProb, addZ, 0, ctr);
    


  // Gate leakage reduction only needed on data after the last CNOT
  if (corLeak == gate) {
    lkgRedCircuit (data4, pProb, puProb, pdProb, NULL);
  }

}

bool Lattice::circSyndromeZZZZ_6 (int xLoc, int yLoc, double pProb, double qProb, double puProb, double pdProb, bool *lkgFlagA, bool *lkgFlagD){

  // Locate the data qubit in the down direction and the ancilla
  Cell *data4, *ancilla;
  if (corLeak == quick) {
    data4   = cells[xLoc][yLoc];
    ancilla = cells[xLoc][(yLoc+1)%ySize];
  } else {
    data4   = cells[xLoc][(yLoc+1)%ySize];
    ancilla = cells[xLoc][yLoc];
  }
  bool result = false;

  // Obtain result by reading the acnilla
  if (ancilla->ErrCurrent->ops[0] == X || ancilla->ErrCurrent->ops[0] == Y)
    result = true;

  // If the ancilla is leaked, the syndrome is also true
  if (!lkgDetect && ancilla->ErrCurrent->ops[0] == L) {
    result = true;
  }

  // If leakage is detected
  *lkgFlagA = (ancilla->ErrCurrent->ops[0] == L);
  if (lkgDetect && ancilla->ErrCurrent->ops[0] == L) {
    result = false;
  }

  // Meanwhile data leaks again...
  data4->ErrCurrent->identityErrLkg (pProb, zero, zero, addZ, 0, ctr);

  // The data needs to be moved back before next round of simulation
  if (corLeak == quick)
    ancilla->ErrCurrent->ops[0] = data4->ErrCurrent->ops[0];

  // Here we do leakage reduction for the circuit
  if (corLeak == circuit) {
    lkgRedCircuit (data4, pProb, puProb, pdProb, lkgFlagD);
  }


  // The syndrome is incorrect with probability p
  if (debug) cout << "Calculating ZZZZ syndrome at xLoc=" << xLoc << " yLoc=" << yLoc << ":"<< endl;
  if (debug) {cout << "  ancilla:"; ancilla->ErrCurrent->printState();}
  if (!lkgDetect || ancilla->ErrCurrent->ops[0] != L) {
    if(ctr->usePath) {
      if (ctr->modeMap[ctr->locCtr] == 1) {
        result = !result;
        if (debug) cout << "syndrome error due to measurement occured" << endl;
      }
    } else if (!ctr->usePath && !ctr->countPath) {
      // The syndrome is incorrect with probability q.
      double r = myRandom();
      if (r < qProb) {
        result = !result;
        if (debug) cout << "  syndrome error due to measurement occured" << endl;
      }
    }
    ctr->locModes.push_back(2);
    ctr->locCauses.push_back(MeasErr);
    ctr->locProbs.push_back(qProb);
    ctr->locCtr++;
  }
  if (debug) cout << "  syndrome=" << result << endl;

  return result;
}


// Call minimum weight perfect matching and correct errors.
void Lattice::callMatching (pauli errType, vector<Vertex*> & vertices, vector<MyEdge*> & edges){

  // Call perfect matching
  struct PerfectMatching::Options options;
  options.verbose = false;
  options.update_duals_before = false;
  options.update_duals_after  = true;
  PerfectMatching *pm = new PerfectMatching((int)vertices.size(),(int)edges.size());
  log << "matching: " << vertices.size() << " vertices, " << edges.size() << " edges" << endl;
  if (debug) cout << "Starting metching with |V|=" << vertices.size() << " and |E|=" << edges.size() << "." << endl; 
  for (int i = 0; i < (int) edges.size(); i++) {
    pm->AddEdge(edges[i]->id1,edges[i]->id2,edges[i]->dist);
  }
  pm->options = options;
  pm->Solve();
  
  // Call the correction subroutine for all matched vertices  
  for (int i = 0; i < (int)vertices.size(); i++ ) {
    // i and j are matched
    int j = pm->GetMatch(i);
    log << "matching: " << i << " and " << j << " are matched" << endl;
    // correct the error unless this is an error on the time axis only 
    if ( i < j && (vertices[i]->x != vertices[j]->x || vertices[i]->y != vertices[j]->y) )
      correctLine (errType, vertices[i]->x, vertices[i]->y, vertices[j]->x, vertices[j]->y);
  }
  delete pm;
}


// Given the locations of two syndromes, correct errors on a line connecting them.
void Lattice::correctLine (pauli errType, int x1Loc, int y1Loc, int x2Loc, int y2Loc) {

  // Initialize x1, y1, x2, y2 so that we can correct errors on a line in 
  // the SE or SW direction coming from coordinates 1 to coordinates 2
  int x1, y1, x2, y2;
  if ( (y1Loc < y2Loc && y2Loc - y1Loc < ySize - y2Loc + y1Loc) || 
       (y2Loc < y1Loc && y1Loc - y2Loc > ySize - y1Loc + y2Loc) ) {
    x1 = x1Loc; y1 = y1Loc; x2 = x2Loc; y2 = y2Loc;
  } else {
    x1 = x2Loc; y1 = y2Loc; x2 = x1Loc; y2 = y1Loc;
  } 

  // Determine if we need to go to SE or SW
  bool east;
  if ( (x1 < x2 && x2 - x1 <  xSize - x2 + x1) ||
       (x2 < x1 && x1 - x2 >  xSize - x1 + x2) )
    east = true;
  else
    east = false;

  if (debug) cout << "Correcting line x1=" << x1 << " y1=" << y1 << " x2=" << x2 << " y2=" << y2 << " east=" << east << ":" << endl;
  log << "correctLine: line from (" << x1 << "," << y1 << ") to (" << x2 << "," << y2 << "), east=" << east << endl;

  // correct in the southbound direction first
  while (y1 != y2) {
    y1 = (y1 + 2) % (ySize);
    correct (errType, x1, (y1-1+ySize) % ySize);
  }

  // if correction in SE direction
  if (east) {
    while (x1 != x2) {
      x1 = (x1 + 2) % (xSize);
      correct (errType, (x1-1+xSize) % xSize, y1);
    }
  // if correction in SW direction
  } else {
    while (x1 != x2) {
      x1 = (x1 - 2 + xSize) % (xSize);
      correct (errType, (x1+1) % xSize, y1);
    }
  }
}


//  Correct a single X or Z error at the specified location.
void Lattice::correct (pauli errType, int xLoc, int yLoc) {

  assert ((errType == X || errType == Z) && xLoc < xSize && yLoc < ySize);

  pauli p;
  if (errType == X)
    p = X;
  else
    p = Z;

  // Correct the error 
  if (debug) cout << " -> correcting " << p << " at x=" << xLoc << " y=" << yLoc << endl;
  log << "correct: " << p << " at x=" << xLoc << " y=" << yLoc << endl;
  cells[xLoc][yLoc]->ErrCurrent->ops[0] = cells[xLoc][yLoc]->ErrCurrent->ops[0] * p;

}


// Generates logical operator X1, X2, Z1, or Z2 on this lattice.
Operator* Lattice::getLogical(string whichOp) {

  // Treat Y1 and Y2 as a product of X and Z
  if (whichOp == "Y1") {
    Operator *myX1 = getLogical ("X1");
    Operator *myZ1 = getLogical ("Z1");
    Operator *myOp = new Operator ((*myX1) * (*myZ1));
    delete myX1; delete myZ1;
    return myOp; 
  } else if (whichOp == "Y2") {
    Operator *myX2 = getLogical ("X2");
    Operator *myZ2 = getLogical ("Z2");
    Operator *myOp = new Operator ((*myX2) * (*myZ2));
    delete myX2; delete myZ2;
    return myOp;
  }

  Operator *myOp = new Operator (xSize*ySize/2);

  for (int i = 0, k = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      if ((i%2 && j%2) || (!(i%2) && !(j%2)))
        continue;
      if (whichOp == "X1") {
        if (j == 1) {
          myOp->ops[k] = X;
        }
      } else if (whichOp == "Z1") {
        if (i == 0) {
          myOp->ops[k] = Z;
        }
      } else if (whichOp == "X2") {
        if (i == 1) {
          myOp->ops[k] = X;
        }
      } else if (whichOp == "Z2") {
        if (j == 0) {
          myOp->ops[k] = Z;
        }
      } else {
        assert (0);
      }
      k++;
    }
  }

  return myOp;
}


// Examimes ErrCurrent to determine if error correction succeeds.
bool Lattice::success(void) {

  // Sanity check: make sure ErrCurrect commutes with all stabilizers
  Operator XXXX(4), ZZZZ(4);
  XXXX.ops[0] = X;
  XXXX.ops[1] = X;
  XXXX.ops[2] = X;
  XXXX.ops[3] = X;
  ZZZZ.ops[0] = Z;
  ZZZZ.ops[1] = Z;
  ZZZZ.ops[2] = Z;
  ZZZZ.ops[3] = Z;
  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      // ZZZZ syndrome
      if ( i%2 && j%2 ) {
        Operator E(0);
        E.pushBack (cells[i][j-1]->ErrCurrent->ops[0]);
	E.pushBack (cells[i-1][j]->ErrCurrent->ops[0]);
	E.pushBack (cells[(i+1)%xSize][j]->ErrCurrent->ops[0]);
	E.pushBack (cells[i][(j+1)%ySize]->ErrCurrent->ops[0]);
        assert (ZZZZ.commute(E));
      }
      // XXXX syndrome
      if ( !(i%2) && !(j%2) ) {
	Operator E(0);
	E.pushBack (cells[i][(j-1+ySize)%ySize]->ErrCurrent->ops[0]);
	E.pushBack (cells[(i-1+xSize)%xSize][j]->ErrCurrent->ops[0]);
	E.pushBack (cells[(i+1)%xSize][j]->ErrCurrent->ops[0]);
	E.pushBack (cells[i][(j+1)%ySize]->ErrCurrent->ops[0]);
	assert (XXXX.commute(E));
      }
    }
  }

  // Obtain ErrCurrent 
  Operator O(0);
  for (int i = 0; i < xSize; i++) {
    for (int j = 0; j < ySize; j++) {
      if ((i%2 && j%2) || (!(i%2) && !(j%2)))
        continue;
      O.pushBack (cells[i][j]->ErrCurrent->ops[0]);
    }
  }

  // Generate logical operators X1, X2, Z1, Z2
  Operator *X1 = getLogical("X1");
  Operator *X2 = getLogical("X2");
  Operator *Z1 = getLogical("Z1");
  Operator *Z2 = getLogical("Z2");

  assert (!X1->commute(*Z1) && !X2->commute(*Z2) && X1->commute(*X2) && X1->commute(*Z2) && X2->commute(*Z1) && Z1->commute(*Z2));

  // Determine result of error correction
  bool result = true;
  if (!O.commute(*X1) || !O.commute(*X2) || !O.commute(*Z1) || !O.commute(*Z2))
    result = false;

  delete X1; delete X2; delete Z1; delete Z2;

  return result;
}


// Print the current state of the lattice. The content of ErrCurrent for all qubits is printed. 
void Lattice::printState(void){

  int xCnt = 0, yCnt = 0, zCnt = 0, lCnt = 0;

  cout << "Lattice: " << xSize << " columns and " << ySize << " rows" << endl;
  for (int j = 0; j < ySize; j++) {
    for (int i = 0; i < xSize; i++) {
      if ((i%2 && j%2) || (!(i%2) && !(j%2)))
        cout << " ";
      else
        cells[i][j]->ErrCurrent->printState();
      switch(cells[i][j]->ErrCurrent->ops[0]) {
      case X: xCnt++; break;
      case Y: yCnt++; break;
      case Z: zCnt++; break;
      case L: lCnt++; break;
      default: break;
      }
    }
    cout << endl;
  }
  cout << "Summary: xCnt=" << xCnt << " yCnt=" << yCnt << " zCnt=" << zCnt << " lCnt=" << lCnt << endl;
}

