#ifndef INTERPOLANT_HH
#define INTERPOLANT_HH

#include "tetrahedron.hh"

class Interpolant
{
public:
  Interpolant(const Tetrahedron &t, double (*ff)(const Vector<3> &),
              Vector<3> (*dff)(const Vector<3> &))
  {
    MPoly<3> linear_basis[4];
    for (int i = 0; i < 4; ++i)
      linear_basis[i] = linear_indicator(t.vertex(i),
                                         t.vertex((i + 1) % 4),
                                         t.vertex((i + 2) % 4),
                                         t.vertex((i + 3) % 4));

    linear_interpolant = 0.0;
    for (int v = 0; v < 4; ++v)
      linear_interpolant += ff(t.vertex(v)) * linear_basis[v];

    cubic_interpolant = linear_interpolant;
    for (int at = 0; at < 4; ++at) {
      Vector<3> vertex = t.vertex(at);
      Vector<3> gradient_want = dff(vertex);
      Vector<3> gradient_have = gradient(linear_interpolant, vertex);
      Vector<3> gradient_difference = gradient_want - gradient_have;
      for (int towards = 0; towards < 4; ++towards) {
        if (at == towards) continue;
        Vector<3> edge = t.edge(at, towards);
        cubic_interpolant +=
          dot_product(gradient_difference, edge)
          * linear_basis[at] * linear_basis[at] * linear_basis[towards];
      }
    }

    quartic_interpolant = cubic_interpolant;
    for (int endpoint1 = 0; endpoint1 < 4; ++endpoint1) {
      for (int endpoint2 = endpoint1 + 1; endpoint2 < 4; ++endpoint2) {
        Vector<3> midpoint = t.edgeMidpoint(endpoint1, endpoint2);
        Vector<3> gradient_want = dff(midpoint);
        Vector<3> gradient_have = gradient(cubic_interpolant, midpoint);
        Vector<3> gradient_difference = gradient_want - gradient_have;
        for (int towards = 0; towards < 4; ++towards) {
          if (endpoint1 == towards || endpoint2 == towards) continue;
          Vector<3> normal = t.edgeNormal(endpoint1, endpoint2, towards);
          quartic_interpolant +=
            dot_product(gradient_difference, normal)
            * 4.0 * linear_basis[endpoint1] * linear_basis[endpoint2]
            * linear_basis[towards]
            * (linear_basis[endpoint1] + linear_basis[endpoint2]
               - linear_basis[towards]);
        }
      }
    }

    quintic_interpolant = quartic_interpolant;
    for (int opposite = 0; opposite < 4; ++opposite) {
      int face1 = (opposite + 1) % 4;
      int face2 = (opposite + 2) % 4;
      int face3 = (opposite + 3) % 4;
      Vector<3> center = t.faceCenter(opposite);
      Vector<3> normal = t.faceNormal(opposite);
      Vector<3> gradient_want = dff(center);
      Vector<3> gradient_have = gradient(quartic_interpolant, center);
      Vector<3> gradient_difference = gradient_want - gradient_have;
      quintic_interpolant +=
        dot_product(gradient_difference, normal)
        * 27.0
        * linear_basis[face1] * linear_basis[face2]
        * linear_basis[face3] * linear_basis[opposite]
        * (1.0 - 3.0 * linear_basis[opposite]);
    }
  }
  const MPoly<3> &linear() const {
    return linear_interpolant;
  }
  const MPoly<3> &cubic() const {
    return cubic_interpolant;
  }
  const MPoly<3> &quartic() const {
    return quartic_interpolant;
  }
  const MPoly<3> &quintic() const {
    return quintic_interpolant;
  }

private:
  MPoly<3> linear_interpolant;
  MPoly<3> cubic_interpolant;
  MPoly<3> quartic_interpolant;
  MPoly<3> quintic_interpolant;

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
};

#endif
