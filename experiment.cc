#include <cstdlib>
#include <cstdio>

#include "mpoly.hh"

#define should(b) should_((b), __LINE__)
inline void should_(bool b, int line)
{
  if (!b) {
    printf("Line %d: Something should be true but isn\'t.\n", line);
    exit(1);
  }
}

int main(int argc, char *argv[])
{
  // Polynomial x
  MPoly<3> x = MPoly<3>::var(0);

  // Vector (3,0,0)
  Vector<3> c(0.0);
  c[0] = 3.0;

  should(x(c) == 3.0);

  return 0;
}
