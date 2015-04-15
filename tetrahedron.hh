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

private:
  Vector<3> vertex_data[4];
};

#endif
