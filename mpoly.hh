#ifndef MPOLY_HH
#define MPOLY_HH

#include <map>

#include "vector.hh"

// Double precision floating point polynomials RR[x_1,...,x_N]

template <int N>
class Monomial
{
private:
  Vector<N> powers;
};

template <int N>
class MPoly
{
private:
  std::map<Monomial<N>,double> terms;
};

#endif
