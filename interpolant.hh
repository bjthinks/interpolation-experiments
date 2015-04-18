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

    MPoly<3> indicator[4];
    indicator[0] = linear_indicator(p, q, r, s);
    indicator[1] = linear_indicator(q, p, r, s);
    indicator[2] = linear_indicator(r, p, q, s);
    indicator[3] = linear_indicator(s, p, q, r);

    linear_poly = 0.0;
    for (int v = 0; v < 4; ++v)
      linear_poly += ff(t.vertex(v)) * indicator[v];

    cubic_poly = linear_poly;
    for (int v_from = 0; v_from < 4; ++v_from) {
      for (int v_to = 0; v_to < 4; ++v_to) {
        if (v_from == v_to)
          continue;
        cubic_poly += dot_product(dff(t.vertex(v_from))
                                  - gradient(linear_poly, t.vertex(v_from)),
                                  t.edge(v_from, v_to))
          * vertexGradient(indicator[v_from], indicator[v_to]);
      }
    }

    quartic_poly = cubic_poly;
    for (int e0 = 0; e0 < 4; ++e0) {
      for (int e1 = e0 + 1; e1 < 4; ++e1) {
        for (int to = 0; to < 4; ++to) {
          if (e0 == to || e1 == to) continue;
          Vector<3> midpoint = t.edgeMidpoint(e0, e1);
          quartic_poly += dot_product(dff(midpoint)
                                    - gradient(cubic_poly, midpoint),
                                    t.edgeNormal(e0, e1, to))
          * edgeGradient(indicator[e0], indicator[e1], indicator[to]);
        }
      }
    }

    quintic_poly = quartic_poly;
    for (int f = 0; f < 4; ++f) {
      quintic_poly +=
        dot_product(dff(t.faceCenter(f))
                    - gradient(quartic_poly, t.faceCenter(f)),
                    t.faceNormal(f))
        * faceGradient(indicator[(f + 1) % 4],
                       indicator[(f + 2) % 4],
                       indicator[(f + 3) % 4],
                       indicator[f]);
    }
  }
  const MPoly<3> &linear() const {
    return linear_poly;
  }
  const MPoly<3> &cubic() const {
    return cubic_poly;
  }
  const MPoly<3> &quartic() const {
    return quartic_poly;
  }
  const MPoly<3> &quintic() const {
    return quintic_poly;
  }

private:
  MPoly<3> linear_poly;
  MPoly<3> cubic_poly;
  MPoly<3> quartic_poly;
  MPoly<3> quintic_poly;

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
  static MPoly<3> edgeGradient(MPoly<3> &e1, MPoly<3> &e2, MPoly<3> &to) {
    return 4.0 * e1 * e2 * to * (e1 + e2 - to);
  }
  static MPoly<3> faceGradient(MPoly<3> &f1, MPoly<3> &f2, MPoly<3> &f3,
                               MPoly<3> &to) {
    return 27.0 * f1 * f2 * f3 * to * (1.0 - 3.0 * to);
  }
};

#endif
