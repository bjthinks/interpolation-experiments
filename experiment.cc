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
  if (fabs(diff) < 256.0 * DBL_EPSILON)
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

  MPoly<3> face_qrs = 27.0 * a*b*c*d * (1.0 - 3.0 * a);
  diag(face_qrs, "face_qrs", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> face_prs = 27.0 * a*b*c*d * (1.0 - 3.0 * b);
  diag(face_prs, "face_prs", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> face_pqs = 27.0 * a*b*c*d * (1.0 - 3.0 * c);
  diag(face_pqs, "face_pqs", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> face_pqr = 27.0 * a*b*c*d * (1.0 - 3.0 * d);
  diag(face_pqr, "face_pqr", p, "p", q, "q", r, "r", s, "s");

  MPoly<3> edge_pq_to_r = 4.0 * a*b*c * (a + b - c)
    + 16.0 / 81.0 * face_pqr
    -  8.0 / 27.0 * face_pqs;
  diag(edge_pq_to_r, "edge_pq_to_r", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> edge_pq_to_s = 4.0 * a*b*d * (a + b - d)
    + 16.0 / 81.0 * face_pqs
    -  8.0 / 27.0 * face_pqr;
  diag(edge_pq_to_s, "edge_pq_to_s", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> edge_pr_to_q = 4.0 * a*b*c * (a + c - b)
    + 16.0 / 81.0 * face_pqr
    -  8.0 / 27.0 * face_prs;
  diag(edge_pr_to_q, "edge_pr_to_q", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> edge_pr_to_s = 4.0 * a*c*d * (a + c - d)
    + 16.0 / 81.0 * face_prs
    -  8.0 / 27.0 * face_pqr;
  diag(edge_pr_to_s, "edge_pr_to_s", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> edge_qr_to_p = 4.0 * a*b*c * (b + c - a)
    + 16.0 / 81.0 * face_pqr
    -  8.0 / 27.0 * face_qrs;
  diag(edge_qr_to_p, "edge_qr_to_p", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> edge_qr_to_s = 4.0 * b*c*d * (b + c - d)
    + 16.0 / 81.0 * face_qrs
    -  8.0 / 27.0 * face_pqr;
  diag(edge_qr_to_s, "edge_qr_to_s", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> edge_ps_to_q = 4.0 * a*b*d * (a + d - b)
    + 16.0 / 81.0 * face_pqs
    -  8.0 / 27.0 * face_prs;
  diag(edge_ps_to_q, "edge_ps_to_q", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> edge_ps_to_r = 4.0 * a*c*d * (a + d - c)
    + 16.0 / 81.0 * face_prs
    -  8.0 / 27.0 * face_pqs;
  diag(edge_ps_to_r, "edge_ps_to_r", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> edge_qs_to_p = 4.0 * a*b*d * (b + d - a)
    + 16.0 / 81.0 * face_pqs
    -  8.0 / 27.0 * face_qrs;
  diag(edge_qs_to_p, "edge_qs_to_p", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> edge_qs_to_r = 4.0 * b*c*d * (b + d - c)
    + 16.0 / 81.0 * face_qrs
    -  8.0 / 27.0 * face_pqs;
  diag(edge_qs_to_r, "edge_qs_to_r", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> edge_rs_to_p = 4.0 * a*c*d * (c + d - a)
    + 16.0 / 81.0 * face_prs
    -  8.0 / 27.0 * face_qrs;
  diag(edge_rs_to_p, "edge_rs_to_p", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> edge_rs_to_q = 4.0 * b*c*d * (c + d - b)
    + 16.0 / 81.0 * face_qrs
    -  8.0 / 27.0 * face_prs;
  diag(edge_rs_to_q, "edge_rs_to_q", p, "p", q, "q", r, "r", s, "s");

  MPoly<3> vert_p_to_q = a*a*b
    + 3.0 / 8.0 * (edge_pq_to_r + edge_pq_to_s)
    - 1.0 / 4.0 * (edge_pr_to_q + edge_ps_to_q)
    + 1.0 / 9.0 * (face_pqr + face_pqs - face_prs);
  diag(vert_p_to_q, "vert_p_to_q", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_p_to_r = a*a*c
    + 3.0 / 8.0 * (edge_pr_to_q + edge_pr_to_s)
    - 1.0 / 4.0 * (edge_pq_to_r + edge_ps_to_r)
    + 1.0 / 9.0 * (face_pqr + face_prs - face_pqs);
  diag(vert_p_to_r, "vert_p_to_r", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_p_to_s = a*a*d
    + 3.0 / 8.0 * (edge_ps_to_q + edge_ps_to_r)
    - 1.0 / 4.0 * (edge_pq_to_s + edge_pr_to_s)
    + 1.0 / 9.0 * (face_pqs + face_prs - face_pqr);
  diag(vert_p_to_s, "vert_p_to_s", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_q_to_p = b*b*a
    + 3.0 / 8.0 * (edge_pq_to_r + edge_pq_to_s)
    - 1.0 / 4.0 * (edge_qr_to_p + edge_qs_to_p)
    + 1.0 / 9.0 * (face_pqr + face_pqs - face_qrs);
  diag(vert_q_to_p, "vert_q_to_p", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_q_to_r = b*b*c
    + 3.0 / 8.0 * (edge_qr_to_p + edge_qr_to_s)
    - 1.0 / 4.0 * (edge_pq_to_r + edge_qs_to_r)
    + 1.0 / 9.0 * (face_pqr + face_qrs - face_pqs);
  diag(vert_q_to_r, "vert_q_to_r", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_q_to_s = b*b*d
    + 3.0 / 8.0 * (edge_qs_to_p + edge_qs_to_r)
    - 1.0 / 4.0 * (edge_pq_to_s + edge_qr_to_s)
    + 1.0 / 9.0 * (face_pqs + face_qrs - face_pqr);
  diag(vert_q_to_s, "vert_q_to_s", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_r_to_p = c*c*a
    + 3.0 / 8.0 * (edge_pr_to_q + edge_pr_to_s)
    - 1.0 / 4.0 * (edge_qr_to_p + edge_rs_to_p)
    + 1.0 / 9.0 * (face_pqr + face_prs - face_qrs);
  diag(vert_r_to_p, "vert_r_to_p", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_r_to_q = c*c*b
    + 3.0 / 8.0 * (edge_qr_to_p + edge_qr_to_s)
    - 1.0 / 4.0 * (edge_pr_to_q + edge_rs_to_q)
    + 1.0 / 9.0 * (face_pqr + face_qrs - face_prs);
  diag(vert_r_to_q, "vert_r_to_q", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_r_to_s = c*c*d
    + 3.0 / 8.0 * (edge_rs_to_p + edge_rs_to_q)
    - 1.0 / 4.0 * (edge_pr_to_s + edge_qr_to_s)
    + 1.0 / 9.0 * (face_prs + face_qrs - face_pqr);
  diag(vert_r_to_s, "vert_r_to_s", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_s_to_p = d*d*a
    + 3.0 / 8.0 * (edge_ps_to_q + edge_ps_to_r)
    - 1.0 / 4.0 * (edge_qs_to_p + edge_rs_to_p)
    + 1.0 / 9.0 * (face_pqs + face_prs - face_qrs);
  diag(vert_s_to_p, "vert_s_to_p", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_s_to_q = d*d*b
    + 3.0 / 8.0 * (edge_qs_to_p + edge_qs_to_r)
    - 1.0 / 4.0 * (edge_ps_to_q + edge_rs_to_q)
    + 1.0 / 9.0 * (face_pqs + face_qrs - face_prs);
  diag(vert_s_to_q, "vert_s_to_q", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_s_to_r = d*d*c
    + 3.0 / 8.0 * (edge_rs_to_p + edge_rs_to_q)
    - 1.0 / 4.0 * (edge_ps_to_r + edge_qs_to_r)
    + 1.0 / 9.0 * (face_prs + face_qrs - face_pqs);
  diag(vert_s_to_r, "vert_s_to_r", p, "p", q, "q", r, "r", s, "s");

  MPoly<3> vert_p = a
    + vert_p_to_q + vert_p_to_r + vert_p_to_s
    - vert_q_to_p - vert_r_to_p - vert_s_to_p
    + 1.0 / 2.0 * (edge_pq_to_r + edge_pq_to_s +
                   edge_pr_to_q + edge_pr_to_s +
                   edge_ps_to_q + edge_ps_to_r)
    - edge_qr_to_p - edge_qs_to_p - edge_rs_to_p
    + 1.0 / 3.0 * (face_pqr + face_pqs + face_prs)
    - face_qrs;
  diag(vert_p, "vert_p", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_q = b
    + vert_q_to_p + vert_q_to_r + vert_q_to_s
    - vert_p_to_q - vert_r_to_q - vert_s_to_q
    + 1.0 / 2.0 * (edge_pq_to_r + edge_pq_to_s +
                   edge_qr_to_p + edge_qr_to_s +
                   edge_qs_to_p + edge_qs_to_r)
    - edge_pr_to_q - edge_ps_to_q - edge_rs_to_q
    + 1.0 / 3.0 * (face_pqr + face_pqs + face_qrs)
    - face_prs;
  diag(vert_q, "vert_q", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_r = c
    + vert_r_to_p + vert_r_to_q + vert_r_to_s
    - vert_p_to_r - vert_q_to_r - vert_s_to_r
    + 1.0 / 2.0 * (edge_pr_to_q + edge_pr_to_s +
                   edge_qr_to_p + edge_qr_to_s +
                   edge_rs_to_p + edge_rs_to_q)
    - edge_pq_to_r - edge_ps_to_r - edge_qs_to_r
    + 1.0 / 3.0 * (face_pqr + face_prs + face_qrs)
    - face_pqs;
  diag(vert_r, "vert_r", p, "p", q, "q", r, "r", s, "s");
  MPoly<3> vert_s = d
    + vert_s_to_p + vert_s_to_q + vert_s_to_r
    - vert_p_to_s - vert_q_to_s - vert_r_to_s
    + 1.0 / 2.0 * (edge_ps_to_q + edge_ps_to_r +
                   edge_qs_to_p + edge_qs_to_r +
                   edge_rs_to_p + edge_rs_to_q)
    - edge_pq_to_s - edge_pr_to_s - edge_qr_to_s
    + 1.0 / 3.0 * (face_pqs + face_prs + face_qrs)
    - face_pqr;
  diag(vert_s, "vert_s", p, "p", q, "q", r, "r", s, "s");

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

  MPoly<3> linear1 =
    p_value * vert_p +
    q_value * vert_q +
    r_value * vert_r +
    s_value * vert_s;
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

  MPoly<3> aab = a * a * b;
  MPoly<3> aac = a * a * c;
  MPoly<3> aad = a * a * d;
  MPoly<3> bba = b * b * a;
  MPoly<3> bbc = b * b * c;
  MPoly<3> bbd = b * b * d;
  MPoly<3> cca = c * c * a;
  MPoly<3> ccb = c * c * b;
  MPoly<3> ccd = c * c * d;
  MPoly<3> dda = d * d * a;
  MPoly<3> ddb = d * d * b;
  MPoly<3> ddc = d * d * c;

  test_zero_values(aab, p, q, r, s);
  test_zero_values(aac, p, q, r, s);
  test_zero_values(aad, p, q, r, s);
  test_zero_values(bba, p, q, r, s);
  test_zero_values(bbc, p, q, r, s);
  test_zero_values(bbd, p, q, r, s);
  test_zero_values(cca, p, q, r, s);
  test_zero_values(ccb, p, q, r, s);
  test_zero_values(ccd, p, q, r, s);
  test_zero_values(dda, p, q, r, s);
  test_zero_values(ddb, p, q, r, s);
  test_zero_values(ddc, p, q, r, s);

  test_cubic_gradients(p, aab, aac, aad, q, r, s);
  test_cubic_gradients(q, bba, bbc, bbd, p, r, s);
  test_cubic_gradients(r, cca, ccb, ccd, p, q, s);
  test_cubic_gradients(s, dda, ddb, ddc, p, q, r);

  MPoly<3> eef = e * e * f;
  MPoly<3> eeg = e * e * g;
  MPoly<3> eeh = e * e * h;
  MPoly<3> ffe = f * f * e;
  MPoly<3> ffg = f * f * g;
  MPoly<3> ffh = f * f * h;
  MPoly<3> gge = g * g * e;
  MPoly<3> ggf = g * g * f;
  MPoly<3> ggh = g * g * h;
  MPoly<3> hhe = h * h * e;
  MPoly<3> hhf = h * h * f;
  MPoly<3> hhg = h * h * g;

  test_zero_values(eef, t, q, r, s);
  test_zero_values(eeg, t, q, r, s);
  test_zero_values(eeh, t, q, r, s);
  test_zero_values(ffe, t, q, r, s);
  test_zero_values(ffg, t, q, r, s);
  test_zero_values(ffh, t, q, r, s);
  test_zero_values(gge, t, q, r, s);
  test_zero_values(ggf, t, q, r, s);
  test_zero_values(ggh, t, q, r, s);
  test_zero_values(hhe, t, q, r, s);
  test_zero_values(hhf, t, q, r, s);
  test_zero_values(hhg, t, q, r, s);

  test_cubic_gradients(t, eef, eeg, eeh, q, r, s);
  test_cubic_gradients(q, ffe, ffg, ffh, t, r, s);
  test_cubic_gradients(r, gge, ggf, ggh, t, q, s);
  test_cubic_gradients(s, hhe, hhf, hhg, t, q, r);

  Vector<3> p_gradient = random_vector<3>();
  Vector<3> q_gradient = random_vector<3>();
  Vector<3> r_gradient = random_vector<3>();
  Vector<3> s_gradient = random_vector<3>();
  Vector<3> t_gradient = random_vector<3>();

  MPoly<3> cubic1 = linear1
    + dot_product(p_gradient, q - p) * vert_p_to_q
    + dot_product(p_gradient, r - p) * vert_p_to_r
    + dot_product(p_gradient, s - p) * vert_p_to_s
    + dot_product(q_gradient, p - q) * vert_q_to_p
    + dot_product(q_gradient, r - q) * vert_q_to_r
    + dot_product(q_gradient, s - q) * vert_q_to_s
    + dot_product(r_gradient, p - r) * vert_r_to_p
    + dot_product(r_gradient, q - r) * vert_r_to_q
    + dot_product(r_gradient, s - r) * vert_r_to_s
    + dot_product(s_gradient, p - s) * vert_s_to_p
    + dot_product(s_gradient, q - s) * vert_s_to_q
    + dot_product(s_gradient, r - s) * vert_s_to_r;

  should(double_equal(cubic1(p), p_value));
  should(double_equal(cubic1(q), q_value));
  should(double_equal(cubic1(r), r_value));
  should(double_equal(cubic1(s), s_value));
  should(vector_equal(gradient(cubic1, p), p_gradient));
  should(vector_equal(gradient(cubic1, q), q_gradient));
  should(vector_equal(gradient(cubic1, r), r_gradient));
  should(vector_equal(gradient(cubic1, s), s_gradient));

  MPoly<3> cubic2 = linear2
    + dot_product(t_gradient - gradient(linear2, t), q - t) * eef
    + dot_product(t_gradient - gradient(linear2, t), r - t) * eeg
    + dot_product(t_gradient - gradient(linear2, t), s - t) * eeh
    + dot_product(q_gradient - gradient(linear2, q), t - q) * ffe
    + dot_product(q_gradient - gradient(linear2, q), r - q) * ffg
    + dot_product(q_gradient - gradient(linear2, q), s - q) * ffh
    + dot_product(r_gradient - gradient(linear2, r), t - r) * gge
    + dot_product(r_gradient - gradient(linear2, r), q - r) * ggf
    + dot_product(r_gradient - gradient(linear2, r), s - r) * ggh
    + dot_product(s_gradient - gradient(linear2, s), t - s) * hhe
    + dot_product(s_gradient - gradient(linear2, s), q - s) * hhf
    + dot_product(s_gradient - gradient(linear2, s), r - s) * hhg;

  should(double_equal(cubic2(t), t_value));
  should(double_equal(cubic2(q), q_value));
  should(double_equal(cubic2(r), r_value));
  should(double_equal(cubic2(s), s_value));
  should(vector_equal(gradient(cubic2, t), t_gradient));
  should(vector_equal(gradient(cubic2, q), q_gradient));
  should(vector_equal(gradient(cubic2, r), r_gradient));
  should(vector_equal(gradient(cubic2, s), s_gradient));

#if 0
  Vector<3> qr = (q + r) / 2.0;
  Vector<3> qs = (q + s) / 2.0;
  Vector<3> rs = (r + s) / 2.0;
  Vector<3> pq = (p + q) / 2.0;
  Vector<3> pr = (p + r) / 2.0;
  Vector<3> ps = (p + s) / 2.0;
  Vector<3> tq = (t + q) / 2.0;
  Vector<3> tr = (t + r) / 2.0;
  Vector<3> ts = (t + s) / 2.0;
#endif

  MPoly<3> quartic1 = cubic1
    + 0.0;

  should(double_equal(quartic1(p), p_value));
  should(double_equal(quartic1(q), q_value));
  should(double_equal(quartic1(r), r_value));
  should(double_equal(quartic1(s), s_value));
  should(vector_equal(gradient(quartic1, p), p_gradient));
  should(vector_equal(gradient(quartic1, q), q_gradient));
  should(vector_equal(gradient(quartic1, r), r_gradient));
  should(vector_equal(gradient(quartic1, s), s_gradient));

  MPoly<3> quartic2 = cubic2
    + 0.0;

  should(double_equal(quartic2(t), t_value));
  should(double_equal(quartic2(q), q_value));
  should(double_equal(quartic2(r), r_value));
  should(double_equal(quartic2(s), s_value));
  should(vector_equal(gradient(quartic2, t), t_gradient));
  should(vector_equal(gradient(quartic2, q), q_gradient));
  should(vector_equal(gradient(quartic2, r), r_gradient));
  should(vector_equal(gradient(quartic2, s), s_gradient));

  //should(vector_equal(gradient(quartic1, qr), gradient(quartic2, qr)));
  //should(vector_equal(gradient(quartic1, qs), gradient(quartic2, qs)));
  //should(vector_equal(gradient(quartic1, rs), gradient(quartic2, rs)));
}
