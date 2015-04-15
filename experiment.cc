#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cfloat>

#include "mpoly.hh"

using namespace std;

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
  if (fabs(diff) < 1e-12)
    return true;
  printf("%f %f %f\n", a, b, a-b);
  return false;
}

template <int N>
bool vector_equal(const Vector<N> &x, const Vector<N> &y) {
  for (int i = 0; i < N; ++i)
    if (!double_equal(x[i], y[i]))
      return false;
  return true;
}

void test_zero_values(const MPoly<3> &f,
                      const Vector<3> &p, const Vector<3> &q,
                      const Vector<3> &r, const Vector<3> &s) {
  should(double_equal(f(p), 0.0));
  should(double_equal(f(q), 0.0));
  should(double_equal(f(r), 0.0));
  should(double_equal(f(s), 0.0));
}

void test_cubic_gradients(const Vector<3> &p, const MPoly<3> &f,
                          const MPoly<3> &g, const MPoly<3> &h,
                          const Vector<3> &q, const Vector<3> &r,
                          const Vector<3> &s) {
  should(double_equal(dot_product(gradient(f, p), q - p), 1.0));
  should(double_equal(dot_product(gradient(f, p), r - p), 0.0));
  should(double_equal(dot_product(gradient(f, p), s - p), 0.0));
  should(double_equal(dot_product(gradient(g, p), q - p), 0.0));
  should(double_equal(dot_product(gradient(g, p), r - p), 1.0));
  should(double_equal(dot_product(gradient(g, p), s - p), 0.0));
  should(double_equal(dot_product(gradient(h, p), q - p), 0.0));
  should(double_equal(dot_product(gradient(h, p), r - p), 0.0));
  should(double_equal(dot_product(gradient(h, p), s - p), 1.0));
}

bool is_zero(double x)
{
  return fabs(x) < 1e-6;
}

void show_value(const MPoly<3> &f,  const string &f_name,
                const Vector<3> &p, const string &p_name)
{
  double value = f(p);
  if (!is_zero(value))
    printf("%s(%s) = %f\n", f_name.c_str(), p_name.c_str(), f(p));
}

void show_dd(const MPoly<3> &f,  const string &f_name,
             const Vector<3> &p, const string &p_name,
             const Vector<3> &q, const string &q_name)
{
  double dd_value = dot_product(gradient(f, p), q - p);
  if (!is_zero(dd_value))
    printf("grad(%s)(%s) dot %s-%s = %f\n", f_name.c_str(), p_name.c_str(),
           q_name.c_str(), p_name.c_str(), dd_value);
}

void diag
(const MPoly<3> &f,  const string &f_name,
 const Vector<3> &p, const string &p_name,
 const Vector<3> &q, const string &q_name,
 const Vector<3> &r, const string &r_name,
 const Vector<3> &s, const string &s_name)
{
  show_value(f, f_name, p, p_name);
  show_value(f, f_name, q, q_name);
  show_value(f, f_name, r, r_name);
  show_value(f, f_name, s, s_name);
  show_dd(f, f_name, p, p_name, q, q_name);
  show_dd(f, f_name, p, p_name, r, r_name);
  show_dd(f, f_name, p, p_name, s, s_name);
  show_dd(f, f_name, q, q_name, p, p_name);
  show_dd(f, f_name, q, q_name, r, r_name);
  show_dd(f, f_name, q, q_name, s, s_name);
  show_dd(f, f_name, r, r_name, p, p_name);
  show_dd(f, f_name, r, r_name, q, q_name);
  show_dd(f, f_name, r, r_name, s, s_name);
  show_dd(f, f_name, s, s_name, p, p_name);
  show_dd(f, f_name, s, s_name, q, q_name);
  show_dd(f, f_name, s, s_name, r, r_name);
  show_dd(f, f_name, (p + q) / 2.0, p_name + q_name, r, r_name);
  show_dd(f, f_name, (p + q) / 2.0, p_name + q_name, s, s_name);
  show_dd(f, f_name, (p + r) / 2.0, p_name + r_name, q, q_name);
  show_dd(f, f_name, (p + r) / 2.0, p_name + r_name, s, s_name);
  show_dd(f, f_name, (p + s) / 2.0, p_name + s_name, q, q_name);
  show_dd(f, f_name, (p + s) / 2.0, p_name + s_name, r, r_name);
  show_dd(f, f_name, (q + r) / 2.0, q_name + r_name, p, p_name);
  show_dd(f, f_name, (q + r) / 2.0, q_name + r_name, s, s_name);
  show_dd(f, f_name, (q + s) / 2.0, q_name + s_name, p, p_name);
  show_dd(f, f_name, (q + s) / 2.0, q_name + s_name, r, r_name);
  show_dd(f, f_name, (r + s) / 2.0, r_name + s_name, p, p_name);
  show_dd(f, f_name, (r + s) / 2.0, r_name + s_name, q, q_name);
  show_dd(f, f_name, (p + q + r) / 3.0, p_name + q_name + r_name, s, s_name);
  show_dd(f, f_name, (p + q + s) / 3.0, p_name + q_name + s_name, r, r_name);
  show_dd(f, f_name, (p + r + s) / 3.0, p_name + r_name + s_name, q, q_name);
  show_dd(f, f_name, (q + r + s) / 3.0, q_name + r_name + s_name, p, p_name);
}

template <int N>
Vector<N> random_on_segment(const Vector<N> &x, const Vector<N> &y) {
  double t = double(rand()) / double(RAND_MAX);
  return t * x + (1.0 - t) * y;
}

MPoly<3> faceGradient(MPoly<3> &f1, MPoly<3> &f2, MPoly<3> &f3,
                      MPoly<3> &opp) {
  return 27.0 * f1 * f2 * f3 * opp * (1.0 - 3.0 * opp);
}

MPoly<3> edgeGradient(MPoly<3> &e1, MPoly<3> &e2, MPoly<3> &to,
                      MPoly<3> &opp) {
  return 4.0 * e1 * e2 * to * (e1 + e2 - to)
    + 16.0 / 81.0 * faceGradient(e1, e2, to, opp)
    -  8.0 / 27.0 * faceGradient(e1, e2, opp, to);
}

MPoly<3> vertexGradient(MPoly<3> &v, MPoly<3> &dir,
                        MPoly<3> &opp1, MPoly<3> &opp2) {
  return v * v * dir
    + 3.0 / 8.0 * (edgeGradient(v, dir, opp1, opp2)
                   + edgeGradient(v, dir, opp2, opp1))
    - 1.0 / 4.0 * (edgeGradient(v, opp1, dir, opp2)
                   + edgeGradient(v, opp2, dir, opp1))
    + 1.0 / 9.0 * (faceGradient(v, dir, opp1, opp2)
                   + faceGradient(v, dir, opp2, opp1)
                   - faceGradient(v, opp1, opp2, dir));
}

MPoly<3> vertexValue(MPoly<3> &v, MPoly<3> &opp1,
                     MPoly<3> &opp2, MPoly<3> &opp3) {
  return v
    + vertexGradient(v, opp1, opp2, opp3)
    + vertexGradient(v, opp2, opp1, opp3)
    + vertexGradient(v, opp3, opp1, opp2)
    - vertexGradient(opp1, v, opp2, opp3)
    - vertexGradient(opp2, v, opp1, opp3)
    - vertexGradient(opp3, v, opp1, opp2)
    + 0.5 * (edgeGradient(v, opp1, opp2, opp3)
             + edgeGradient(v, opp1, opp3, opp2)
             + edgeGradient(v, opp2, opp1, opp3)
             + edgeGradient(v, opp2, opp3, opp1)
             + edgeGradient(v, opp3, opp1, opp2)
             + edgeGradient(v, opp3, opp2, opp1))
    - edgeGradient(opp1, opp2, v, opp3)
    - edgeGradient(opp1, opp3, v, opp2)
    - edgeGradient(opp2, opp3, v, opp1)
    + 1.0 / 3.0 * (faceGradient(v, opp1, opp2, opp3)
                   + faceGradient(v, opp1, opp3, opp2)
                   + faceGradient(v, opp2, opp3, opp1))
    - faceGradient(opp1, opp2, opp3, v);
}

int main(int argc, char *argv[]) {
  srand(345987);

  mpoly_tests();
  mpoly_diff_tests();

  // Make two random tetrahedra with a shared face: (p,q,r,s) & (t,q,r,s)

  Vector<3> p = random_vector<3>();
  Vector<3> q = random_vector<3>();
  Vector<3> r = random_vector<3>();
  Vector<3> s = random_vector<3>();
  Vector<3> t = random_vector<3>();

  // Define some linear functions on the tetrahedra

  MPoly<3> a = linear_indicator(p, q, r, s);
  MPoly<3> b = linear_indicator(q, p, r, s);
  MPoly<3> c = linear_indicator(r, p, q, s);
  MPoly<3> d = linear_indicator(s, p, q, r);

  MPoly<3> e = linear_indicator(t, q, r, s);
  MPoly<3> f = linear_indicator(q, t, r, s);
  MPoly<3> g = linear_indicator(r, t, q, s);
  MPoly<3> h = linear_indicator(s, t, q, r);

  // Make sure they have the desired properties

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

  // Next, define a basis of functions for the first tetrahedron,
  // and verify its properties

  MPoly<3> vert_a = vertexValue(a, b, c, d);
  MPoly<3> vert_b = vertexValue(b, a, c, d);
  MPoly<3> vert_c = vertexValue(c, a, b, d);
  MPoly<3> vert_d = vertexValue(d, a, b, c);

  diag(vert_a, "vert_a", p, "p", q, "q", r, "r", s, "s");
  diag(vert_b, "vert_b", p, "p", q, "q", r, "r", s, "s");
  diag(vert_c, "vert_c", p, "p", q, "q", r, "r", s, "s");
  diag(vert_d, "vert_d", p, "p", q, "q", r, "r", s, "s");

  MPoly<3> vert_a_to_b = vertexGradient(a, b, c, d);
  MPoly<3> vert_a_to_c = vertexGradient(a, c, b, d);
  MPoly<3> vert_a_to_d = vertexGradient(a, d, b, c);
  MPoly<3> vert_b_to_a = vertexGradient(b, a, c, d);
  MPoly<3> vert_b_to_c = vertexGradient(b, c, a, d);
  MPoly<3> vert_b_to_d = vertexGradient(b, d, a, c);
  MPoly<3> vert_c_to_a = vertexGradient(c, a, b, d);
  MPoly<3> vert_c_to_b = vertexGradient(c, b, a, d);
  MPoly<3> vert_c_to_d = vertexGradient(c, d, a, b);
  MPoly<3> vert_d_to_a = vertexGradient(d, a, b, c);
  MPoly<3> vert_d_to_b = vertexGradient(d, b, a, c);
  MPoly<3> vert_d_to_c = vertexGradient(d, c, a, b);

  diag(vert_a_to_b, "vert_a_to_b", p, "p", q, "q", r, "r", s, "s");
  diag(vert_a_to_c, "vert_a_to_c", p, "p", q, "q", r, "r", s, "s");
  diag(vert_a_to_d, "vert_a_to_d", p, "p", q, "q", r, "r", s, "s");
  diag(vert_b_to_a, "vert_b_to_a", p, "p", q, "q", r, "r", s, "s");
  diag(vert_b_to_c, "vert_b_to_c", p, "p", q, "q", r, "r", s, "s");
  diag(vert_b_to_d, "vert_b_to_d", p, "p", q, "q", r, "r", s, "s");
  diag(vert_c_to_a, "vert_c_to_a", p, "p", q, "q", r, "r", s, "s");
  diag(vert_c_to_b, "vert_c_to_b", p, "p", q, "q", r, "r", s, "s");
  diag(vert_c_to_d, "vert_c_to_d", p, "p", q, "q", r, "r", s, "s");
  diag(vert_d_to_a, "vert_d_to_a", p, "p", q, "q", r, "r", s, "s");
  diag(vert_d_to_b, "vert_d_to_b", p, "p", q, "q", r, "r", s, "s");
  diag(vert_d_to_c, "vert_d_to_c", p, "p", q, "q", r, "r", s, "s");

  MPoly<3> edge_ab_to_c = edgeGradient(a, b, c, d);
  MPoly<3> edge_ab_to_d = edgeGradient(a, b, d, c);
  MPoly<3> edge_ac_to_b = edgeGradient(a, c, b, d);
  MPoly<3> edge_ac_to_d = edgeGradient(a, c, d, b);
  MPoly<3> edge_bc_to_a = edgeGradient(b, c, a, d);
  MPoly<3> edge_bc_to_d = edgeGradient(b, c, d, a);
  MPoly<3> edge_ad_to_b = edgeGradient(a, d, b, c);
  MPoly<3> edge_ad_to_c = edgeGradient(a, d, c, b);
  MPoly<3> edge_bd_to_a = edgeGradient(b, d, a, c);
  MPoly<3> edge_bd_to_c = edgeGradient(b, d, c, a);
  MPoly<3> edge_cd_to_a = edgeGradient(c, d, a, b);
  MPoly<3> edge_cd_to_b = edgeGradient(c, d, b, a);

  diag(edge_ab_to_c, "edge_ab_to_c", p, "p", q, "q", r, "r", s, "s");
  diag(edge_ab_to_d, "edge_ab_to_d", p, "p", q, "q", r, "r", s, "s");
  diag(edge_ac_to_b, "edge_ac_to_b", p, "p", q, "q", r, "r", s, "s");
  diag(edge_ac_to_d, "edge_ac_to_d", p, "p", q, "q", r, "r", s, "s");
  diag(edge_bc_to_a, "edge_bc_to_a", p, "p", q, "q", r, "r", s, "s");
  diag(edge_bc_to_d, "edge_bc_to_d", p, "p", q, "q", r, "r", s, "s");
  diag(edge_ad_to_b, "edge_ad_to_b", p, "p", q, "q", r, "r", s, "s");
  diag(edge_ad_to_c, "edge_ad_to_c", p, "p", q, "q", r, "r", s, "s");
  diag(edge_bd_to_a, "edge_bd_to_a", p, "p", q, "q", r, "r", s, "s");
  diag(edge_bd_to_c, "edge_bd_to_c", p, "p", q, "q", r, "r", s, "s");
  diag(edge_cd_to_a, "edge_cd_to_a", p, "p", q, "q", r, "r", s, "s");
  diag(edge_cd_to_b, "edge_cd_to_b", p, "p", q, "q", r, "r", s, "s");

  MPoly<3> face_bcd = faceGradient(b, c, d, a);
  MPoly<3> face_acd = faceGradient(a, c, d, b);
  MPoly<3> face_abd = faceGradient(a, b, d, c);
  MPoly<3> face_abc = faceGradient(a, b, c, d);

  diag(face_bcd, "face_bcd", p, "p", q, "q", r, "r", s, "s");
  diag(face_acd, "face_acd", p, "p", q, "q", r, "r", s, "s");
  diag(face_abd, "face_abd", p, "p", q, "q", r, "r", s, "s");
  diag(face_abc, "face_abc", p, "p", q, "q", r, "r", s, "s");

  // Now, do the same for the second tetrahedron

  MPoly<3> vert_e = vertexValue(e, f, g, h);
  MPoly<3> vert_f = vertexValue(f, e, g, h);
  MPoly<3> vert_g = vertexValue(g, e, f, h);
  MPoly<3> vert_h = vertexValue(h, e, f, g);

  diag(vert_e, "vert_e", t, "t", q, "q", r, "r", s, "s");
  diag(vert_f, "vert_f", t, "t", q, "q", r, "r", s, "s");
  diag(vert_g, "vert_g", t, "t", q, "q", r, "r", s, "s");
  diag(vert_h, "vert_h", t, "t", q, "q", r, "r", s, "s");

  MPoly<3> vert_e_to_f = vertexGradient(e, f, g, h);
  MPoly<3> vert_e_to_g = vertexGradient(e, g, f, h);
  MPoly<3> vert_e_to_h = vertexGradient(e, h, f, g);
  MPoly<3> vert_f_to_e = vertexGradient(f, e, g, h);
  MPoly<3> vert_f_to_g = vertexGradient(f, g, e, h);
  MPoly<3> vert_f_to_h = vertexGradient(f, h, e, g);
  MPoly<3> vert_g_to_e = vertexGradient(g, e, f, h);
  MPoly<3> vert_g_to_f = vertexGradient(g, f, e, h);
  MPoly<3> vert_g_to_h = vertexGradient(g, h, e, f);
  MPoly<3> vert_h_to_e = vertexGradient(h, e, f, g);
  MPoly<3> vert_h_to_f = vertexGradient(h, f, e, g);
  MPoly<3> vert_h_to_g = vertexGradient(h, g, e, f);

  diag(vert_e_to_f, "vert_e_to_f", t, "t", q, "q", r, "r", s, "s");
  diag(vert_e_to_g, "vert_e_to_g", t, "t", q, "q", r, "r", s, "s");
  diag(vert_e_to_h, "vert_e_to_h", t, "t", q, "q", r, "r", s, "s");
  diag(vert_f_to_e, "vert_f_to_e", t, "t", q, "q", r, "r", s, "s");
  diag(vert_f_to_g, "vert_f_to_g", t, "t", q, "q", r, "r", s, "s");
  diag(vert_f_to_h, "vert_f_to_h", t, "t", q, "q", r, "r", s, "s");
  diag(vert_g_to_e, "vert_g_to_e", t, "t", q, "q", r, "r", s, "s");
  diag(vert_g_to_f, "vert_g_to_f", t, "t", q, "q", r, "r", s, "s");
  diag(vert_g_to_h, "vert_g_to_h", t, "t", q, "q", r, "r", s, "s");
  diag(vert_h_to_e, "vert_h_to_e", t, "t", q, "q", r, "r", s, "s");
  diag(vert_h_to_f, "vert_h_to_f", t, "t", q, "q", r, "r", s, "s");
  diag(vert_h_to_g, "vert_h_to_g", t, "t", q, "q", r, "r", s, "s");

  MPoly<3> edge_ef_to_g = edgeGradient(e, f, g, h);
  MPoly<3> edge_ef_to_h = edgeGradient(e, f, h, g);
  MPoly<3> edge_eg_to_f = edgeGradient(e, g, f, h);
  MPoly<3> edge_eg_to_h = edgeGradient(e, g, h, f);
  MPoly<3> edge_fg_to_e = edgeGradient(f, g, e, h);
  MPoly<3> edge_fg_to_h = edgeGradient(f, g, h, e);
  MPoly<3> edge_eh_to_f = edgeGradient(e, h, f, g);
  MPoly<3> edge_eh_to_g = edgeGradient(e, h, g, f);
  MPoly<3> edge_fh_to_e = edgeGradient(f, h, e, g);
  MPoly<3> edge_fh_to_g = edgeGradient(f, h, g, e);
  MPoly<3> edge_gh_to_e = edgeGradient(g, h, e, f);
  MPoly<3> edge_gh_to_f = edgeGradient(g, h, f, e);

  diag(edge_ef_to_g, "edge_ef_to_g", t, "t", q, "q", r, "r", s, "s");
  diag(edge_ef_to_h, "edge_ef_to_h", t, "t", q, "q", r, "r", s, "s");
  diag(edge_eg_to_f, "edge_eg_to_f", t, "t", q, "q", r, "r", s, "s");
  diag(edge_eg_to_h, "edge_eg_to_h", t, "t", q, "q", r, "r", s, "s");
  diag(edge_fg_to_e, "edge_fg_to_e", t, "t", q, "q", r, "r", s, "s");
  diag(edge_fg_to_h, "edge_fg_to_h", t, "t", q, "q", r, "r", s, "s");
  diag(edge_eh_to_f, "edge_eh_to_f", t, "t", q, "q", r, "r", s, "s");
  diag(edge_eh_to_g, "edge_eh_to_g", t, "t", q, "q", r, "r", s, "s");
  diag(edge_fh_to_e, "edge_fh_to_e", t, "t", q, "q", r, "r", s, "s");
  diag(edge_fh_to_g, "edge_fh_to_g", t, "t", q, "q", r, "r", s, "s");
  diag(edge_gh_to_e, "edge_gh_to_e", t, "t", q, "q", r, "r", s, "s");
  diag(edge_gh_to_f, "edge_gh_to_f", t, "t", q, "q", r, "r", s, "s");

  MPoly<3> face_fgh = faceGradient(f, g, h, e);
  MPoly<3> face_egh = faceGradient(e, g, h, f);
  MPoly<3> face_efh = faceGradient(e, f, h, g);
  MPoly<3> face_efg = faceGradient(e, f, g, h);

  diag(face_fgh, "face_fgh", t, "t", q, "q", r, "r", s, "s");
  diag(face_egh, "face_egh", t, "t", q, "q", r, "r", s, "s");
  diag(face_efh, "face_efh", t, "t", q, "q", r, "r", s, "s");
  diag(face_efg, "face_efg", t, "t", q, "q", r, "r", s, "s");

  // Make up some values for the interpolant to have at the vertices

  double p_value = random_double();
  double q_value = random_double();
  double r_value = random_double();
  double s_value = random_double();
  double t_value = random_double();

  // Make interpolants that have these values

  MPoly<3> values1 =
    p_value * vert_a +
    q_value * vert_b +
    r_value * vert_c +
    s_value * vert_d;
  MPoly<3> values2 =
    t_value * vert_e +
    q_value * vert_f +
    r_value * vert_g +
    s_value * vert_h;

  // And check that they are correct

  should(double_equal(values1(p), p_value));
  should(double_equal(values1(q), q_value));
  should(double_equal(values1(r), r_value));
  should(double_equal(values1(s), s_value));
  should(double_equal(values2(t), t_value));
  should(double_equal(values2(q), q_value));
  should(double_equal(values2(r), r_value));
  should(double_equal(values2(s), s_value));

  // Now make up some gradients at each vertex

  Vector<3> p_gradient = random_vector<3>();
  Vector<3> q_gradient = random_vector<3>();
  Vector<3> r_gradient = random_vector<3>();
  Vector<3> s_gradient = random_vector<3>();
  Vector<3> t_gradient = random_vector<3>();

  // And make interpolants that have these vertex gradients

  MPoly<3> vgrads1 = values1
    + dot_product(p_gradient, q - p) * vert_a_to_b
    + dot_product(p_gradient, r - p) * vert_a_to_c
    + dot_product(p_gradient, s - p) * vert_a_to_d
    + dot_product(q_gradient, p - q) * vert_b_to_a
    + dot_product(q_gradient, r - q) * vert_b_to_c
    + dot_product(q_gradient, s - q) * vert_b_to_d
    + dot_product(r_gradient, p - r) * vert_c_to_a
    + dot_product(r_gradient, q - r) * vert_c_to_b
    + dot_product(r_gradient, s - r) * vert_c_to_d
    + dot_product(s_gradient, p - s) * vert_d_to_a
    + dot_product(s_gradient, q - s) * vert_d_to_b
    + dot_product(s_gradient, r - s) * vert_d_to_c;

  MPoly<3> vgrads2 = values2
    + dot_product(t_gradient, q - t) * vert_e_to_f
    + dot_product(t_gradient, r - t) * vert_e_to_g
    + dot_product(t_gradient, s - t) * vert_e_to_h
    + dot_product(q_gradient, t - q) * vert_f_to_e
    + dot_product(q_gradient, r - q) * vert_f_to_g
    + dot_product(q_gradient, s - q) * vert_f_to_h
    + dot_product(r_gradient, t - r) * vert_g_to_e
    + dot_product(r_gradient, q - r) * vert_g_to_f
    + dot_product(r_gradient, s - r) * vert_g_to_h
    + dot_product(s_gradient, t - s) * vert_h_to_e
    + dot_product(s_gradient, q - s) * vert_h_to_f
    + dot_product(s_gradient, r - s) * vert_h_to_g;

  // And check that they are correct

  should(double_equal(vgrads1(p), p_value));
  should(double_equal(vgrads1(q), q_value));
  should(double_equal(vgrads1(r), r_value));
  should(double_equal(vgrads1(s), s_value));
  should(vector_equal(gradient(vgrads1, p), p_gradient));
  should(vector_equal(gradient(vgrads1, q), q_gradient));
  should(vector_equal(gradient(vgrads1, r), r_gradient));
  should(vector_equal(gradient(vgrads1, s), s_gradient));

  should(double_equal(vgrads2(t), t_value));
  should(double_equal(vgrads2(q), q_value));
  should(double_equal(vgrads2(r), r_value));
  should(double_equal(vgrads2(s), s_value));
  should(vector_equal(gradient(vgrads2, t), t_gradient));
  should(vector_equal(gradient(vgrads2, q), q_gradient));
  should(vector_equal(gradient(vgrads2, r), r_gradient));
  should(vector_equal(gradient(vgrads2, s), s_gradient));

  // Now define the edge midpoints, and make up gradients to approximate
  // at each of them

  // Shared edges

  Vector<3> qr = (q + r) / 2.0;
  //Vector<3> qs = (q + s) / 2.0;
  //Vector<3> rs = (r + s) / 2.0;

  Vector<3> qr_gradient = random_vector<3>();
  //Vector<3> qs_gradient = random_vector<3>();
  //Vector<3> rs_gradient = random_vector<3>();
#if 0
  // Edges unique to first tetrahedron

  Vector<3> pq = (p + q) / 2.0;
  Vector<3> pr = (p + r) / 2.0;
  Vector<3> ps = (p + s) / 2.0;

  Vector<3> pq_gradient = random_vector<3>();
  Vector<3> pr_gradient = random_vector<3>();
  Vector<3> ps_gradient = random_vector<3>();

  // Edges unique to second tetrahedron

  Vector<3> tq = (t + q) / 2.0;
  Vector<3> tr = (t + r) / 2.0;
  Vector<3> ts = (t + s) / 2.0;

  Vector<3> tq_gradient = random_vector<3>();
  Vector<3> tr_gradient = random_vector<3>();
  Vector<3> ts_gradient = random_vector<3>();
#endif
  // And make interpolants that have these vertex gradients

  MPoly<3> egrads1 = vgrads1
    + dot_product(qr_gradient, p - qr) * edge_bc_to_a
    + dot_product(qr_gradient, s - qr) * edge_bc_to_d;

  should(double_equal(egrads1(p), p_value));
  should(double_equal(egrads1(q), q_value));
  should(double_equal(egrads1(r), r_value));
  should(double_equal(egrads1(s), s_value));
  should(vector_equal(gradient(egrads1, p), p_gradient));
  should(vector_equal(gradient(egrads1, q), q_gradient));
  should(vector_equal(gradient(egrads1, r), r_gradient));
  should(vector_equal(gradient(egrads1, s), s_gradient));

  MPoly<3> egrads2 = vgrads2
    + dot_product(qr_gradient, t - qr) * edge_fg_to_e
    + dot_product(qr_gradient, s - qr) * edge_fg_to_h;

  should(double_equal(egrads2(t), t_value));
  should(double_equal(egrads2(q), q_value));
  should(double_equal(egrads2(r), r_value));
  should(double_equal(egrads2(s), s_value));
  should(vector_equal(gradient(egrads2, t), t_gradient));
  should(vector_equal(gradient(egrads2, q), q_gradient));
  should(vector_equal(gradient(egrads2, r), r_gradient));
  should(vector_equal(gradient(egrads2, s), s_gradient));

  should(vector_equal(gradient(egrads1, qr), gradient(egrads2, qr)));
  //should(vector_equal(gradient(egrads1, qs), gradient(egrads2, qs)));
  //should(vector_equal(gradient(egrads1, rs), gradient(egrads2, rs)));
}
