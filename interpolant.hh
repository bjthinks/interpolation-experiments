#ifndef INTERPOLANT_HH
#define INTERPOLANT_HH

#include "tetrahedron.hh"

class Interpolant
{
public:
  Interpolant(const Tetrahedron &t_)
    : t(t_)
  {
  }

private:
  const Tetrahedron &t;
};

#endif
