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
    for (int i = 0; i < N; ++i) {
      if (powers[i] < rhs.powers[i])
        return true;
      if (powers[i] > rhs.powers[i])
        return false;
    }
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
class MPoly : public VectorSpace<double, MPoly<N> >
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
    for (c_iter i = terms.begin(); i != terms.end(); ++i)
      value += (i->second) * (i->first)(p);
    return value;
  }

  MPoly<N> operator-() const
  {
    MPoly<N> f;
    for (c_iter i = terms.begin(); i != terms.end(); ++i)
      f.terms[(i->first)] = -(i->second);
    return f;
  }

  const MPoly<N> &operator+=(const MPoly &rhs)
  {
    for (c_iter i = rhs.terms.begin(); i != rhs.terms.end(); ++i) {
      iter j = terms.find(i->first);
      if (j == terms.end())
        terms[i->first] = i->second;
      else {
        j->second += i->second;
        if (j->second == 0.0)
          terms.erase(i->first);
      }
    }
    return *this;
  }

  const MPoly<N> &operator*=(double rhs)
  {
    if (rhs == 0.0) {
      terms = std::map<Monomial<N>,double>();
      return *this;
    }

    for (iter i = terms.begin(); i != terms.end(); ++i) {
      i->second *= rhs;
    }
    return *this;
  }

private:
  std::map<Monomial<N>,double> terms;
  typedef typename std::map<Monomial<N>,double>::iterator iter;
  typedef typename std::map<Monomial<N>,double>::const_iterator c_iter;
};

#endif
