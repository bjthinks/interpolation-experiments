#ifndef INTERPOLANT_HH
#define INTERPOLANT_HH

#include "tetrahedron.hh"

class Interpolant
{
public:
  Interpolant(const Tetrahedron &t, double (*ff)(const Vector<3> &),
              Vector<3> (*dff)(const Vector<3> &))
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

    cubic_poly = linear_poly
      + dot_product(dff(p) - gradient(linear_poly, p),
                    t.edge(0, 1)) * vertexGradient(a, b)
      + dot_product(dff(p) - gradient(linear_poly, p),
                    t.edge(0, 2)) * vertexGradient(a, c)
      + dot_product(dff(p) - gradient(linear_poly, p),
                    t.edge(0, 3)) * vertexGradient(a, d)
      + dot_product(dff(q) - gradient(linear_poly, q),
                    t.edge(1, 0)) * vertexGradient(b, a)
      + dot_product(dff(q) - gradient(linear_poly, q),
                    t.edge(1, 2)) * vertexGradient(b, c)
      + dot_product(dff(q) - gradient(linear_poly, q),
                    t.edge(1, 3)) * vertexGradient(b, d)
      + dot_product(dff(r) - gradient(linear_poly, r),
                    t.edge(2, 0)) * vertexGradient(c, a)
      + dot_product(dff(r) - gradient(linear_poly, r),
                    t.edge(2, 1)) * vertexGradient(c, b)
      + dot_product(dff(r) - gradient(linear_poly, r),
                    t.edge(2, 3)) * vertexGradient(c, d)
      + dot_product(dff(s) - gradient(linear_poly, s),
                    t.edge(3, 0)) * vertexGradient(d, a)
      + dot_product(dff(s) - gradient(linear_poly, s),
                    t.edge(3, 1)) * vertexGradient(d, b)
      + dot_product(dff(s) - gradient(linear_poly, s),
                    t.edge(3, 2)) * vertexGradient(d, c);
  }
  const MPoly<3> &linear() const {
    return linear_poly;
  }
  const MPoly<3> &cubic() const {
    return cubic_poly;
  }

private:
  MPoly<3> linear_poly;
  MPoly<3> cubic_poly;

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
  static Vector<3> gradient(const MPoly<3> &f, const Vector<3> &p) {
    Vector<3> g;
    for (int i = 0; i < 3; ++i)
      g[i] = f.diff(i)(p);
    return g;
  }
  static MPoly<3> vertexGradient(MPoly<3> &v, MPoly<3> &to) {
    return v * v * to;
  }
};

#endif
