#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cfloat>

#include "mpoly.hh"

#define should(b) should_((b), __LINE__)
inline void should_(bool b, int line) {
  if (!b) {
    printf("Line %d: Something should be true but isn\'t.\n", line);
    exit(1);
  }
}

void mpoly_tests() {
  // Variables
  MPoly<3> x = MPoly<3>::var(0);
  MPoly<3> y = MPoly<3>::var(1);
  MPoly<3> z = MPoly<3>::var(2);

  // Vectors a = (3,0,0), b = (0,5,0), c = (0,0,7)
  Vector<3> a(0.0), b(0.0), c(0.0);
  a[0] = 3.0;
  b[1] = 5.0;
  c[2] = 7.0;

  should(x(a) == 3.0);
  should(x(b) == 0.0);
  should(x(c) == 0.0);
  should(y(a) == 0.0);
  should(y(b) == 5.0);
  should(y(c) == 0.0);
  should(z(a) == 0.0);
  should(z(b) == 0.0);
  should(z(c) == 7.0);

  // Polynomial negation
  should((-x)(a) == -3.0);
  should((-y)(b) == -5.0);
  should((-z)(c) == -7.0);

  // Polynomial +=
  MPoly<3> x_plus_y(x);
  x_plus_y += y;
  should(x_plus_y(a) == 3.0);
  should(x_plus_y(b) == 5.0);
  should(x_plus_y(c) == 0.0);

  // Polynomial *= constant
  MPoly<3> x_times_4(x);
  x_times_4 *= 4.0;
  should(x_times_4(a) == 12.0);
  should(x_times_4(b) == 0.0);
  should(x_times_4(c) == 0.0);

  // Now MPoly<N> should be a VectorSpace<double,MPoly<N> >
  should((+x)(a) == 3.0);
  should((+x)(b) == 0.0);
  should((+x)(c) == 0.0);
  should((x + y + z)(a) == 3.0);
  should((x + y + z)(b) == 5.0);
  should((x + y + z)(c) == 7.0);
  should((x + y - z)(a) == 3.0);
  should((x + y - z)(b) == 5.0);
  should((x + y - z)(c) == -7.0);
  should((2.0 * x - 3.0 * y + 4.0 * z)(a) == 6.0);
  should((2.0 * x - 3.0 * y + 4.0 * z)(b) == -15.0);
  should((2.0 * x - 3.0 * y + 4.0 * z)(c) == 28.0);

  // Explicit construction from doubles
  MPoly<3> one(1.0);
  should(one(a) == 1.0);
  should(one(b) == 1.0);
  should(one(c) == 1.0);

  // Now MPoly<N> should be an Algebra<double,MPoly<N> >
  should((x + 1.0)(a) == 4.0);
  should((y + 1.0)(a) == 1.0);
  should((z + 1.0)(a) == 1.0);
  MPoly<3> v;
  v = 2.0;
  should(v(a) == 2.0);
  should(v(b) == 2.0);
  should(v(c) == 2.0);

  MPoly<3> xyz;
  xyz = x * y * z;
  should(xyz(a) == 0.0);
  should(xyz(b) == 0.0);
  should(xyz(c) == 0.0);
  should(xyz(a+b) == 0.0);
  should(xyz(b+c) == 0.0);
  should(xyz(a+c) == 0.0);
  should(xyz(a + b + c) == 105.0);

  // And a CommutativeAlgebra<double,MPoly<N> >
  MPoly<3> product(1.0);
  product *= x;
  product *= y;
  product *= z;
  should(product(a) == 0.0);
  should(product(b) == 0.0);
  should(product(c) == 0.0);
  should(product(a+b) == 0.0);
  should(product(b+c) == 0.0);
  should(product(a+c) == 0.0);
  should(product(a + b + c) == 105.0);
}

template <int N>
Vector<N> gradient(const MPoly<N> &f, const Vector<N> &p) {
  Vector<N> g;
  for (int i = 0; i < N; ++i)
    g[i] = f.diff(i)(p);
  return g;
}

void mpoly_diff_tests() {
  MPoly<1> t = MPoly<1>::var(0);
  MPoly<1> k(1.0);
  should(k.diff(0)(Vector<1>(0.0)) == 0.0);
  should(k.diff(0)(Vector<1>(3.0)) == 0.0);
  should(t.diff(0)(Vector<1>(0.0)) == 1.0);
  should(t.diff(0)(Vector<1>(3.0)) == 1.0);
  should((t * t).diff(0)(Vector<1>(0.0)) == 0.0);
  should((t * t).diff(0)(Vector<1>(3.0)) == 6.0);

  // Variables
  MPoly<3> x = MPoly<3>::var(0);
  MPoly<3> y = MPoly<3>::var(1);
  MPoly<3> z = MPoly<3>::var(2);

  // Test point
  Vector<3> p;
  p[0] = p[1] = p[2] = 1.0;

  MPoly<3> f = x * x * y + 3.0 * y * y * z - 5.0 * x * z * z + 7.0 * x * x * x;
  should(f(p) == 6.0);
  Vector<3> result;
  result[0] = 18.0;
  result[1] = 7.0;
  result[2] = -7.0;
  should(f.diff(0)(p) == result[0]);
  should(f.diff(1)(p) == result[1]);
  should(f.diff(2)(p) == result[2]);
  should(gradient(f, p) == result);
}

double random_double() {
  return 2.0 * double(rand()) / double(RAND_MAX) - 1.0;
}

template <int N>
Vector<N> random_vector() {
  Vector<N> foo;
  for (int i = 0; i < N; ++i)
    foo[i] = random_double();
  return foo;
}

MPoly<3> linear_indicator(const Vector<3> &one,
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

bool double_equal(double a, double b) {
  double diff = a - b;
  if (fabs(diff) < 32.0 * DBL_EPSILON)
    return true;
  return false;
}

int main(int argc, char *argv[]) {
  srand(345987);

  mpoly_tests();

  MPoly<3> x(MPoly<3>::var(0));
  MPoly<3> y(MPoly<3>::var(1));
  MPoly<3> z(MPoly<3>::var(2));
  MPoly<3> k(1.0);

  Vector<3> p = random_vector<3>();
  Vector<3> q = random_vector<3>();
  Vector<3> r = random_vector<3>();
  Vector<3> s = random_vector<3>();
  Vector<3> t = random_vector<3>();

  MPoly<3> a = linear_indicator(p, q, r, s);
  MPoly<3> b = linear_indicator(q, p, r, s);
  MPoly<3> c = linear_indicator(r, p, q, s);
  MPoly<3> d = linear_indicator(s, p, q, r);

  should(double_equal(a(p), 1.0));
  should(double_equal(a(q), 0.0));
  should(double_equal(a(r), 0.0));
  should(double_equal(a(s), 0.0));
  should(double_equal(b(p), 0.0));
  should(double_equal(b(q), 1.0));
  should(double_equal(b(r), 0.0));
  should(double_equal(b(s), 0.0));
  should(double_equal(c(p), 0.0));
  should(double_equal(c(q), 0.0));
  should(double_equal(c(r), 1.0));
  should(double_equal(c(s), 0.0));
  should(double_equal(d(p), 0.0));
  should(double_equal(d(q), 0.0));
  should(double_equal(d(r), 0.0));
  should(double_equal(d(s), 1.0));

  MPoly<3> e = linear_indicator(t, q, r, s);
  MPoly<3> f = linear_indicator(q, t, r, s);
  MPoly<3> g = linear_indicator(r, t, q, s);
  MPoly<3> h = linear_indicator(s, t, q, r);

  should(double_equal(e(t), 1.0));
  should(double_equal(e(q), 0.0));
  should(double_equal(e(r), 0.0));
  should(double_equal(e(s), 0.0));
  should(double_equal(f(t), 0.0));
  should(double_equal(f(q), 1.0));
  should(double_equal(f(r), 0.0));
  should(double_equal(f(s), 0.0));
  should(double_equal(g(t), 0.0));
  should(double_equal(g(q), 0.0));
  should(double_equal(g(r), 1.0));
  should(double_equal(g(s), 0.0));
  should(double_equal(h(t), 0.0));
  should(double_equal(h(q), 0.0));
  should(double_equal(h(r), 0.0));
  should(double_equal(h(s), 1.0));

  double p_value = random_double();
  double q_value = random_double();
  double r_value = random_double();
  double s_value = random_double();
  double t_value = random_double();

  MPoly<3> linear1 = p_value * a + q_value * b + r_value * c + s_value * d;
  MPoly<3> linear2 = t_value * e + q_value * f + r_value * g + s_value * h;

  should(double_equal(linear1(p), p_value));
  should(double_equal(linear1(q), q_value));
  should(double_equal(linear1(r), r_value));
  should(double_equal(linear1(s), s_value));
  should(double_equal(linear2(t), t_value));
  should(double_equal(linear2(q), q_value));
  should(double_equal(linear2(r), r_value));
  should(double_equal(linear2(s), s_value));

  mpoly_diff_tests();
}
