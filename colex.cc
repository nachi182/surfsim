/*! \file colex.cc
 *  \brief enumerate combinations
 *  \author Andrew W Cross <awcross@us.ibm.com>
 *  \date June 2005
 */

#include "colex.hh"

/*! constructor
 * \param ni total number of objects
 * \param ki number to choose
 */
colex::colex(int ni, int ki)
{
  n = ni;
  k = ki;	
  v = new int [n+2];
  for(i = 0; i < n+2; i++) v[i] = 0;
  r = new int [k+1];
  r[k] = 0;
}

/*! destructor
 */
colex::~colex()
{
  delete [] v;
  delete [] r;
}

/*! From Constructive Combinatorics by Stanton and White.
 * \return array (length k) while combinations, else NULL
 */
int *colex::next(void)
{
  if(v[1] == 0)
  {
    for(i = 1; i <= k; i++)
    {
      v[i] = i;
      r[i-1] = v[i];
    }
    v[k+1] = n + 1;
    return r;
  }
  else if(v[1] < n - k + 1)
  {
    j = 0;
    do { j++; } while (v[j+1] <= v[j] + 1);
    v[j]++;
    for(i = 1; i < j; i++) v[i] = i;
    for(i = 1; i <= k; i++) r[i-1] = v[i];
    return r;
  }
  else
    return 0;
}

/*! return the number of combinations
 * \return number of combinations
 */
unsigned long colex::size(void)
{
  return nchoosek(n,k);
}

/*! compute the number of combinations n choose k
 * \param n number of objects
 * \param k number to choose
 * \return n choose k
 * \bug this should be done with big-int package
 */
unsigned long nchoosek(int n, int k)
{
  unsigned long c = 1UL;
  if(2*k > n) k = n - k;
  for(unsigned long i = 1 ; (int)i <= k; i++)
  {
    unsigned long f = n - (i - 1);
    if( f % i == 0 ) { f /= i; assert( c*f >= c ); c *= f; }
    else if( c % i == 0 ) { c /= i; assert( c*f >= c ); c *= f; }
    else { assert( c*f >= c ); c *= f; c /= i; }
  }
  return c;
}
