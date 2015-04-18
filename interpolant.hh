#ifndef INTERPOLANT_HH
#define INTERPOLANT_HH

#include "tetrahedron.hh"

class Interpolant
{
public:
  Interpolant(const Tetrahedron &t, double (*ff)(const Vector<3> &))
  {
    Vector<3> p = t.vertex(0);
    Vector<3> q = t.vertex(1);
    Vector<3> r = t.vertex(2);
    Vector<3> s = t.vertex(3);

    MPoly<3> a = linear_indicator(p, q, r, s);
    MPoly<3> b = linear_indicator(q, p, r, s);
    MPoly<3> c = linear_indicator(r, p, q, s);
    MPoly<3> d = linear_indicator(s, p, q, r);

    linear_poly = ff(p) * a + ff(q) * b + ff(r) * c + ff(s) * d;
  }
  const MPoly<3> &linear() const {
    return linear_poly;
  }

private:
  MPoly<3> linear_poly;

  static MPoly<3> linear_indicator(const Vector<3> &one,
                                   const Vector<3> &zero1,
                                   const Vector<3> &zero2,
                                   const Vector<3> &zero3) {
    Vector<3> zero_normal = cross_product(zero2 - zero1, zero3 - zero1);
    MPoly<3> pre_result
      = zero_normal[0] * MPoly<3>::var(0)
      + zero_normal[1] * MPoly<3>::var(1)
      + zero_normal[2] * MPoly<3>::var(2);
    return (pre_result - pre_result(zero1))
      / (pre_result(one) - pre_result(zero1));
  }
};

#endif
