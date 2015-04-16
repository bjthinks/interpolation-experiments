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
    return vertex_data[v];
  }
  Vector<3> edge(int from, int to) const
  {
    return vertex(to) - vertex(from);
  }
  Vector<3> faceCenter(int f) const
  {
    int p = (f + 1) % 4;
    int q = (f + 2) % 4;
    int r = (f + 3) % 4;
    return (vertex(p) + vertex(q) + vertex(r)) / 3.0;
  }
  Vector<3> faceNormalUnscaled(int f) const
  {
    int p = (f + 1) % 4;
    int q = (f + 2) % 4;
    int r = (f + 3) % 4;
    return cross_product(edge(r, p), edge(r, q));
  }
  Vector<3> faceNormal(int f) const
  {
    return projection(vertex(f) - faceCenter(f), faceNormalUnscaled(f));
  }

private:
  Vector<3> vertex_data[4];

  static Vector<3> projection(const Vector<3> &vec, const Vector<3> &onto)
  {
    return dot_product(vec, onto) / dot_product(onto, onto) * onto;
  }
};

#endif
