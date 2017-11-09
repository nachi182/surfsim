/*! \file colex.hh
 *  \brief enumerate combinations
 *  \author Andrew W Cross <awcross@us.ibm.com>
 *  \date 13 June 2005
 */

#ifndef __COLEX__
#define __COLEX__

#include <cassert>

using namespace std;

unsigned long nchoosek(int n, int k);

/*! class for enumerating combinations
 */
class colex
{
  protected:
    int i, j, k, n; // for looping, storing n, k
    int *v; // for building combination 
    int *r; // for returning combination

  public:
    colex(int ni, int ki);
    ~colex();
    int *next(void);
    unsigned long size(void);
};

#endif
