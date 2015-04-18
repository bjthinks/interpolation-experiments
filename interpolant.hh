#ifndef INTERPOLANT_HH
#define INTERPOLANT_HH

#include "tetrahedron.hh"

class Interpolant
{
public:
  Interpolant(const Tetrahedron &t_)
    : t(t_)
  {
    Vector<3> p = t.vertex(0);
    Vector<3> q = t.vertex(0);
    Vector<3> r = t.vertex(0);
    Vector<3> s = t.vertex(0);

    MPoly<3> a = linear_indicator(p, q, r, s);
    MPoly<3> b = linear_indicator(q, p, r, s);
    MPoly<3> c = linear_indicator(r, p, q, s);
    MPoly<3> d = linear_indicator(s, p, q, r);
  }

private:
  const Tetrahedron &t;
};

#endif
