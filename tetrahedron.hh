#ifndef TETRAHEDRON_HH
#define TETRAHEDRON_HH

class Tetrahedron
{
public:
  Tetrahedron(const Vector<3> &vertex0, const Vector<3> &vertex1,
              const Vector<3> &vertex2, const Vector<3> &vertex3)
  {
    vertex_data[0] = vertex0;
    vertex_data[1] = vertex1;
    vertex_data[2] = vertex2;
    vertex_data[3] = vertex3;
  }
  const Vector<3> &vertex(int v) const
  {
    if (v < 0 || v >= 4)
      throw 0;
    return vertex_data[v];
  }
  Vector<3> faceNormal(int v) const
  {
    int p = (v + 1) % 4;
    int q = (v + 2) % 4;
    int r = (v + 3) % 4;
    return project(vertex(v) - (vertex(p) + vertex(q) + vertex(r)) / 3.0,
                   cross_product(vertex(p) - vertex(r), vertex(q) - vertex(r)));
  }

private:
  Vector<3> vertex_data[4];
  static Vector<3> project(const Vector<3> &vec, const Vector<3> &onto) {
    return dot_product(vec, onto) / dot_product(onto, onto) * onto;
  }
};

#endif
