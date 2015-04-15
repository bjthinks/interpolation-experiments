#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cfloat>

#include "mpoly.hh"
#include "tetrahedron.hh"

using namespace std;

#define should(b) should_((b), __LINE__)
inline void should_(bool b, int line) {
  if (!b) {
    printf("Line %d: Something should be true but isn\'t.\n", line);
    exit(1);
  }
}

static void mpoly_tests() {
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
static Vector<N> gradient(const MPoly<N> &f, const Vector<N> &p) {
  Vector<N> g;
  for (int i = 0; i < N; ++i)
    g[i] = f.diff(i)(p);
  return g;
}

static void mpoly_diff_tests() {
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

static double random_double() {
  return 2.0 * double(rand()) / double(RAND_MAX) - 1.0;
}

template <int N>
static Vector<N> random_vector() {
  Vector<N> foo;
  for (int i = 0; i < N; ++i)
    foo[i] = random_double();
  return foo;
}

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

static bool double_equal(double a, double b) {
  double diff = a - b;
  if (fabs(diff) < 1e-12)
    return true;
  printf("%f %f %f\n", a, b, a-b);
  return false;
}

template <int N>
static bool vector_equal(const Vector<N> &x, const Vector<N> &y) {
  for (int i = 0; i < N; ++i)
    if (!double_equal(x[i], y[i]))
      return false;
  return true;
}

static Vector<3> project(const Vector<3> &vec, const Vector<3> &onto) {
  return dot_product(vec, onto) / dot_product(onto, onto) * onto;
}

static Vector<3> perp(const Vector<3> &vec, const Vector<3> &away) {
  return vec - project(vec, away);
}

template <int N>
static Vector<N> random_on_segment(const Vector<N> &x, const Vector<N> &y) {
  double t = double(rand()) / double(RAND_MAX);
  return t * x + (1.0 - t) * y;
}

static MPoly<3> faceGradient(MPoly<3> &f1, MPoly<3> &f2, MPoly<3> &f3,
                             MPoly<3> &to) {
  return 27.0 * f1 * f2 * f3 * to * (1.0 - 3.0 * to);
}

static MPoly<3> edgeGradient(MPoly<3> &e1, MPoly<3> &e2, MPoly<3> &to) {
  return 4.0 * e1 * e2 * to * (e1 + e2 - to);
}

static MPoly<3> vertexGradient(MPoly<3> &v, MPoly<3> &to) {
  return v * v * to;
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

  // Make up some values for the interpolant to have at the vertices

  double p_value = random_double();
  double q_value = random_double();
  double r_value = random_double();
  double s_value = random_double();
  double t_value = random_double();

  // Make interpolants that have these values

  MPoly<3> values1 = p_value * a + q_value * b + r_value * c + s_value * d;
  MPoly<3> values2 = t_value * e + q_value * f + r_value * g + s_value * h;

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
    + dot_product(p_gradient - gradient(values1, p),
                  q - p) * vertexGradient(a, b)
    + dot_product(p_gradient - gradient(values1, p),
                  r - p) * vertexGradient(a, c)
    + dot_product(p_gradient - gradient(values1, p),
                  s - p) * vertexGradient(a, d)
    + dot_product(q_gradient - gradient(values1, q),
                  p - q) * vertexGradient(b, a)
    + dot_product(q_gradient - gradient(values1, q),
                  r - q) * vertexGradient(b, c)
    + dot_product(q_gradient - gradient(values1, q),
                  s - q) * vertexGradient(b, d)
    + dot_product(r_gradient - gradient(values1, r),
                  p - r) * vertexGradient(c, a)
    + dot_product(r_gradient - gradient(values1, r),
                  q - r) * vertexGradient(c, b)
    + dot_product(r_gradient - gradient(values1, r),
                  s - r) * vertexGradient(c, d)
    + dot_product(s_gradient - gradient(values1, s),
                  p - s) * vertexGradient(d, a)
    + dot_product(s_gradient - gradient(values1, s),
                  q - s) * vertexGradient(d, b)
    + dot_product(s_gradient - gradient(values1, s),
                  r - s) * vertexGradient(d, c);

  MPoly<3> vgrads2 = values2
    + dot_product(t_gradient - gradient(values2, t),
                  q - t) * vertexGradient(e, f)
    + dot_product(t_gradient - gradient(values2, t),
                  r - t) * vertexGradient(e, g)
    + dot_product(t_gradient - gradient(values2, t),
                  s - t) * vertexGradient(e, h)
    + dot_product(q_gradient - gradient(values2, q),
                  t - q) * vertexGradient(f, e)
    + dot_product(q_gradient - gradient(values2, q),
                  r - q) * vertexGradient(f, g)
    + dot_product(q_gradient - gradient(values2, q),
                  s - q) * vertexGradient(f, h)
    + dot_product(r_gradient - gradient(values2, r),
                  t - r) * vertexGradient(g, e)
    + dot_product(r_gradient - gradient(values2, r),
                  q - r) * vertexGradient(g, f)
    + dot_product(r_gradient - gradient(values2, r),
                  s - r) * vertexGradient(g, h)
    + dot_product(s_gradient - gradient(values2, s),
                  t - s) * vertexGradient(h, e)
    + dot_product(s_gradient - gradient(values2, s),
                  q - s) * vertexGradient(h, f)
    + dot_product(s_gradient - gradient(values2, s),
                  r - s) * vertexGradient(h, g);

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
  Vector<3> qs = (q + s) / 2.0;
  Vector<3> rs = (r + s) / 2.0;

  Vector<3> qr_gradient = random_vector<3>();
  Vector<3> qs_gradient = random_vector<3>();
  Vector<3> rs_gradient = random_vector<3>();

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

  // Make interpolants that approximate these edge gradients

  MPoly<3> egrads1 = vgrads1
    + dot_product(pq_gradient - gradient(vgrads1, pq),
                  perp(r - pq, p - q)) * edgeGradient(a, b, c)
    + dot_product(pq_gradient - gradient(vgrads1, pq),
                  perp(s - pq, p - q)) * edgeGradient(a, b, d)
    + dot_product(pr_gradient - gradient(vgrads1, pr),
                  perp(q - pr, p - r)) * edgeGradient(a, c, b)
    + dot_product(pr_gradient - gradient(vgrads1, pr),
                  perp(s - pr, p - r)) * edgeGradient(a, c, d)
    + dot_product(ps_gradient - gradient(vgrads1, ps),
                  perp(q - ps, p - s)) * edgeGradient(a, d, b)
    + dot_product(ps_gradient - gradient(vgrads1, ps),
                  perp(r - ps, p - s)) * edgeGradient(a, d, c)
    + dot_product(qr_gradient - gradient(vgrads1, qr),
                  perp(p - qr, q - r)) * edgeGradient(b, c, a)
    + dot_product(qr_gradient - gradient(vgrads1, qr),
                  perp(s - qr, q - r)) * edgeGradient(b, c, d)
    + dot_product(qs_gradient - gradient(vgrads1, qs),
                  perp(p - qs, q - s)) * edgeGradient(b, d, a)
    + dot_product(qs_gradient - gradient(vgrads1, qs),
                  perp(r - qs, q - s)) * edgeGradient(b, d, c)
    + dot_product(rs_gradient - gradient(vgrads1, rs),
                  perp(p - rs, r - s)) * edgeGradient(c, d, a)
    + dot_product(rs_gradient - gradient(vgrads1, rs),
                  perp(q - rs, r - s)) * edgeGradient(c, d, b);

  MPoly<3> egrads2 = vgrads2
    + dot_product(tq_gradient - gradient(vgrads2, tq),
                  perp(r - tq, t - q)) * edgeGradient(e, f, g)
    + dot_product(tq_gradient - gradient(vgrads2, tq),
                  perp(s - tq, t - q)) * edgeGradient(e, f, h)
    + dot_product(tr_gradient - gradient(vgrads2, tr),
                  perp(q - tr, t - r)) * edgeGradient(e, g, f)
    + dot_product(tr_gradient - gradient(vgrads2, tr),
                  perp(s - tr, t - r)) * edgeGradient(e, g, h)
    + dot_product(ts_gradient - gradient(vgrads2, ts),
                  perp(q - ts, t - s)) * edgeGradient(e, h, f)
    + dot_product(ts_gradient - gradient(vgrads2, ts),
                  perp(r - ts, t - s)) * edgeGradient(e, h, g)
    + dot_product(qr_gradient - gradient(vgrads2, qr),
                  perp(t - qr, q - r)) * edgeGradient(f, g, e)
    + dot_product(qr_gradient - gradient(vgrads2, qr),
                  perp(s - qr, q - r)) * edgeGradient(f, g, h)
    + dot_product(qs_gradient - gradient(vgrads2, qs),
                  perp(t - qs, q - s)) * edgeGradient(f, h, e)
    + dot_product(qs_gradient - gradient(vgrads2, qs),
                  perp(r - qs, q - s)) * edgeGradient(f, h, g)
    + dot_product(rs_gradient - gradient(vgrads2, rs),
                  perp(t - rs, r - s)) * edgeGradient(g, h, e)
    + dot_product(rs_gradient - gradient(vgrads2, rs),
                  perp(q - rs, r - s)) * edgeGradient(g, h, f);

  // And check that they are correct

  should(double_equal(egrads1(p), p_value));
  should(double_equal(egrads1(q), q_value));
  should(double_equal(egrads1(r), r_value));
  should(double_equal(egrads1(s), s_value));
  should(vector_equal(gradient(egrads1, p), p_gradient));
  should(vector_equal(gradient(egrads1, q), q_gradient));
  should(vector_equal(gradient(egrads1, r), r_gradient));
  should(vector_equal(gradient(egrads1, s), s_gradient));
  should(vector_equal(perp(pq_gradient, p - q),
                      perp(gradient(egrads1, pq), p - q)));
  should(vector_equal(perp(pr_gradient, p - r),
                      perp(gradient(egrads1, pr), p - r)));
  should(vector_equal(perp(ps_gradient, p - s),
                      perp(gradient(egrads1, ps), p - s)));
  should(vector_equal(perp(qr_gradient, q - r),
                      perp(gradient(egrads1, qr), q - r)));
  should(vector_equal(perp(qs_gradient, q - s),
                      perp(gradient(egrads1, qs), q - s)));
  should(vector_equal(perp(rs_gradient, r - s),
                      perp(gradient(egrads1, rs), r - s)));

  should(double_equal(egrads2(t), t_value));
  should(double_equal(egrads2(q), q_value));
  should(double_equal(egrads2(r), r_value));
  should(double_equal(egrads2(s), s_value));
  should(vector_equal(gradient(egrads2, t), t_gradient));
  should(vector_equal(gradient(egrads2, q), q_gradient));
  should(vector_equal(gradient(egrads2, r), r_gradient));
  should(vector_equal(gradient(egrads2, s), s_gradient));
  should(vector_equal(perp(tq_gradient, t - q),
                      perp(gradient(egrads2, tq), t - q)));
  should(vector_equal(perp(tr_gradient, t - r),
                      perp(gradient(egrads2, tr), t - r)));
  should(vector_equal(perp(ts_gradient, t - s),
                      perp(gradient(egrads2, ts), t - s)));
  should(vector_equal(perp(qr_gradient, q - r),
                      perp(gradient(egrads2, qr), q - r)));
  should(vector_equal(perp(qs_gradient, q - s),
                      perp(gradient(egrads2, qs), q - s)));
  should(vector_equal(perp(rs_gradient, r - s),
                      perp(gradient(egrads2, rs), r - s)));

  // In particular, edge gradients of the two interpolants
  // should match exactly.

  should(vector_equal(gradient(egrads1, qr), gradient(egrads2, qr)));
  should(vector_equal(gradient(egrads1, qs), gradient(egrads2, qs)));
  should(vector_equal(gradient(egrads1, rs), gradient(egrads2, rs)));

  // Now define the face centers, and make up gradients to approximate
  // at each of them

  // Shared face

  Vector<3> qrs = (q + r + s) / 3.0;
  Vector<3> qrs_gradient = random_vector<3>();

  // Faces unique to first tetrahedron

  Vector<3> pqr = (p + q + r) / 3.0;
  Vector<3> pqs = (p + q + s) / 3.0;
  Vector<3> prs = (p + r + s) / 3.0;

  Vector<3> pqr_gradient = random_vector<3>();
  Vector<3> pqs_gradient = random_vector<3>();
  Vector<3> prs_gradient = random_vector<3>();

  // Faces unique to second tetrahedron

  Vector<3> tqr = (t + q + r) / 3.0;
  Vector<3> tqs = (t + q + s) / 3.0;
  Vector<3> trs = (t + r + s) / 3.0;

  Vector<3> tqr_gradient = random_vector<3>();
  Vector<3> tqs_gradient = random_vector<3>();
  Vector<3> trs_gradient = random_vector<3>();

  // Make interpolants that approximate these face gradients

  MPoly<3> fgrads1 = egrads1
    + dot_product(pqr_gradient - gradient(egrads1, pqr),
                  project(s - pqr, cross_product(p - r, q - r))
                  ) * faceGradient(a, b, c, d)
    + dot_product(pqs_gradient - gradient(egrads1, pqs),
                  project(r - pqs, cross_product(p - s, q - s))
                  ) * faceGradient(a, b, d, c)
    + dot_product(prs_gradient - gradient(egrads1, prs),
                  project(q - prs, cross_product(p - s, r - s))
                  ) * faceGradient(a, c, d, b)
    + dot_product(qrs_gradient - gradient(egrads1, qrs),
                  project(p - qrs, cross_product(q - s, r - s))
                  ) * faceGradient(b, c, d, a);

  MPoly<3> fgrads2 = egrads2
    + dot_product(tqr_gradient - gradient(egrads2, tqr),
                  project(s - tqr, cross_product(t - r, q - r))
                  ) * faceGradient(e, f, g, h)
    + dot_product(tqs_gradient - gradient(egrads2, tqs),
                  project(r - tqs, cross_product(t - s, q - s))
                  ) * faceGradient(e, f, h, g)
    + dot_product(trs_gradient - gradient(egrads2, trs),
                  project(q - trs, cross_product(t - s, r - s))
                  ) * faceGradient(e, g, h, f)
    + dot_product(qrs_gradient - gradient(egrads2, qrs),
                  project(t - qrs, cross_product(q - s, r - s))
                  ) * faceGradient(f, g, h, e);

  // And check that they are correct

  should(double_equal(fgrads1(p), p_value));
  should(double_equal(fgrads1(q), q_value));
  should(double_equal(fgrads1(r), r_value));
  should(double_equal(fgrads1(s), s_value));
  should(vector_equal(gradient(fgrads1, p), p_gradient));
  should(vector_equal(gradient(fgrads1, q), q_gradient));
  should(vector_equal(gradient(fgrads1, r), r_gradient));
  should(vector_equal(gradient(fgrads1, s), s_gradient));
  should(vector_equal(perp(pq_gradient, p - q),
                      perp(gradient(fgrads1, pq), p - q)));
  should(vector_equal(perp(pr_gradient, p - r),
                      perp(gradient(fgrads1, pr), p - r)));
  should(vector_equal(perp(ps_gradient, p - s),
                      perp(gradient(fgrads1, ps), p - s)));
  should(vector_equal(perp(qr_gradient, q - r),
                      perp(gradient(fgrads1, qr), q - r)));
  should(vector_equal(perp(qs_gradient, q - s),
                      perp(gradient(fgrads1, qs), q - s)));
  should(vector_equal(perp(rs_gradient, r - s),
                      perp(gradient(fgrads1, rs), r - s)));
  should(vector_equal(project(pqr_gradient, cross_product(p - r, q - r)),
                      project(gradient(fgrads1, pqr),
                              cross_product(p - r, q - r))));
  should(vector_equal(project(pqs_gradient, cross_product(p - s, q - s)),
                      project(gradient(fgrads1, pqs),
                              cross_product(p - s, q - s))));
  should(vector_equal(project(prs_gradient, cross_product(p - s, r - s)),
                      project(gradient(fgrads1, prs),
                              cross_product(p - s, r - s))));
  should(vector_equal(project(qrs_gradient, cross_product(q - s, r - s)),
                      project(gradient(fgrads1, qrs),
                              cross_product(q - s, r - s))));

  should(double_equal(fgrads2(t), t_value));
  should(double_equal(fgrads2(q), q_value));
  should(double_equal(fgrads2(r), r_value));
  should(double_equal(fgrads2(s), s_value));
  should(vector_equal(gradient(fgrads2, t), t_gradient));
  should(vector_equal(gradient(fgrads2, q), q_gradient));
  should(vector_equal(gradient(fgrads2, r), r_gradient));
  should(vector_equal(gradient(fgrads2, s), s_gradient));
  should(vector_equal(perp(tq_gradient, t - q),
                      perp(gradient(fgrads2, tq), t - q)));
  should(vector_equal(perp(tr_gradient, t - r),
                      perp(gradient(fgrads2, tr), t - r)));
  should(vector_equal(perp(ts_gradient, t - s),
                      perp(gradient(fgrads2, ts), t - s)));
  should(vector_equal(perp(qr_gradient, q - r),
                      perp(gradient(fgrads2, qr), q - r)));
  should(vector_equal(perp(qs_gradient, q - s),
                      perp(gradient(fgrads2, qs), q - s)));
  should(vector_equal(perp(rs_gradient, r - s),
                      perp(gradient(fgrads2, rs), r - s)));
  should(vector_equal(project(tqr_gradient, cross_product(t - r, q - r)),
                      project(gradient(fgrads2, tqr),
                              cross_product(t - r, q - r))));
  should(vector_equal(project(tqs_gradient, cross_product(t - s, q - s)),
                      project(gradient(fgrads2, tqs),
                              cross_product(t - s, q - s))));
  should(vector_equal(project(trs_gradient, cross_product(t - s, r - s)),
                      project(gradient(fgrads2, trs),
                              cross_product(t - s, r - s))));
  should(vector_equal(project(qrs_gradient, cross_product(q - s, r - s)),
                      project(gradient(fgrads2, qrs),
                              cross_product(q - s, r - s))));

  // Edge gradients of the two interpolants should still match exactly,

  should(vector_equal(gradient(fgrads1, qr), gradient(fgrads2, qr)));
  should(vector_equal(gradient(fgrads1, qs), gradient(fgrads2, qs)));
  should(vector_equal(gradient(fgrads1, rs), gradient(fgrads2, rs)));

  // And so should the face gradient

  should(vector_equal(gradient(fgrads1, qrs), gradient(fgrads2, qrs)));
}
