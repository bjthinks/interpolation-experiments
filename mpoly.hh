#ifndef MPOLY_HH
#define MPOLY_HH

#include <map>

#include "vector.hh"

// Double precision floating point polynomials RR[x_1,...,x_N]

template <int N>
class Monomial
{
public:
  Monomial()
  {
    for (int i = 0; i < N; ++i)
      powers[i] = 0;
  }

  static Monomial<N> var(int n)
  {
    Monomial<N> x;
    x.powers[n] = 1;
    return x;
  }

  bool operator<(const Monomial<N> &rhs) const
  {
    for (int i = 0; i < N; ++i)
      if (powers[i] < rhs.powers[i])
        return true;
    return false;
  }

  double operator()(const Vector<N> &p) const
  {
    double value = 1.0;
    for (int i = 0; i < N; ++i)
      value *= pow(p[i], powers[i]);
    return value;
  }

private:
  int powers[N];
};

template <int N>
class MPoly
{
public:
  static MPoly<N> var(int n)
  {
    Monomial<N> x = Monomial<N>::var(n);
    MPoly<N> f;
    f.terms[x] = 1.0;
    return f;
  }

  double operator()(const Vector<N> &p) const
  {
    double value = 0.0;
    typename std::map<Monomial<N>,double>::const_iterator i;
    for (i = terms.begin(); i != terms.end(); ++i)
      value += (i->second) * (i->first)(p);
    return value;
  }

private:
  std::map<Monomial<N>,double> terms;
};

#endif
