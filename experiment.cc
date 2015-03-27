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

void show_value(const MPoly<3> &f, const char *f_name,
                const Vector<3> &p, const char *p_name)
{
  printf("%s(%s) = %f\n", f_name, p_name, f(p));
}

void show_dd(const MPoly<3> &f, const char *f_name,
             const Vector<3> &p, const char *p_name,
             const Vector<3> &q, const char *q_name)
{
  printf("grad(%s)(%s) dot %s-%s = %f\n", f_name, p_name, q_name, p_name,
         dot_product(gradient(f, p), q - p));
}

void diag
(const MPoly<3> &f, const char *f_name,
 const Vector<3> &p, const char *p_name, const Vector<3> &q, const char *q_name,
 const Vector<3> &r, const char *r_name, const Vector<3> &s, const char *s_name)
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
#if 0
  printf("aabc\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(aabc, qr), p - qr),
         dot_product(gradient(aabc, qr), s - qr),
         dot_product(gradient(aabc, qs), p - qs),
         dot_product(gradient(aabc, qs), r - qs),
         dot_product(gradient(aabc, rs), p - rs),
         dot_product(gradient(aabc, rs), q - rs));
#endif
}

template <int N>
Vector<N> random_on_segment(const Vector<N> &x, const Vector<N> &y) {
  double t = double(rand()) / double(RAND_MAX);
  return t * x + (1.0 - t) * y;
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

  diag(a, "a", p, "p", q, "q", r, "r", s, "s");

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
    + dot_product(p_gradient - gradient(linear1, p), q - p) * aab
    + dot_product(p_gradient - gradient(linear1, p), r - p) * aac
    + dot_product(p_gradient - gradient(linear1, p), s - p) * aad
    + dot_product(q_gradient - gradient(linear1, q), p - q) * bba
    + dot_product(q_gradient - gradient(linear1, q), r - q) * bbc
    + dot_product(q_gradient - gradient(linear1, q), s - q) * bbd
    + dot_product(r_gradient - gradient(linear1, r), p - r) * cca
    + dot_product(r_gradient - gradient(linear1, r), q - r) * ccb
    + dot_product(r_gradient - gradient(linear1, r), s - r) * ccd
    + dot_product(s_gradient - gradient(linear1, s), p - s) * dda
    + dot_product(s_gradient - gradient(linear1, s), q - s) * ddb
    + dot_product(s_gradient - gradient(linear1, s), r - s) * ddc;

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

  Vector<3> qr = (q + r) / 2.0;
  Vector<3> qs = (q + s) / 2.0;
  Vector<3> rs = (r + s) / 2.0;
  Vector<3> pq = (p + q) / 2.0;
  Vector<3> pr = (p + r) / 2.0;
  Vector<3> ps = (p + s) / 2.0;
  Vector<3> tq = (t + q) / 2.0;
  Vector<3> tr = (t + r) / 2.0;
  Vector<3> ts = (t + s) / 2.0;

  MPoly<3> aabc = a * a * b * c;
  MPoly<3> aabd = a * a * b * d;
  MPoly<3> aacd = a * a * c * d;
  MPoly<3> bbac = b * b * a * c;
  MPoly<3> bbad = b * b * a * d;
  MPoly<3> bbcd = b * b * c * d;
  MPoly<3> ccab = c * c * a * b;
  MPoly<3> ccad = c * c * a * d;
  MPoly<3> ccbd = c * c * b * d;
  MPoly<3> ddab = d * d * a * b;
  MPoly<3> ddac = d * d * a * c;
  MPoly<3> ddbc = d * d * b * c;

  MPoly<3> t1_qr_p = 4.0 * (- aabc + bbac + ccab);
  MPoly<3> t1_pr_q = 4.0 * (+ aabc - bbac + ccab);
  MPoly<3> t1_pq_r = 4.0 * (+ aabc + bbac - ccab);
  MPoly<3> t1_qs_p = 4.0 * (- aabd + bbad + ddab);
  MPoly<3> t1_ps_q = 4.0 * (+ aabd - bbad + ddab);
  MPoly<3> t1_pq_s = 4.0 * (+ aabd + bbad - ddab);
  MPoly<3> t1_rs_p = 4.0 * (- aacd + ccad + ddac);
  MPoly<3> t1_ps_r = 4.0 * (+ aacd - ccad + ddac);
  MPoly<3> t1_pr_s = 4.0 * (+ aacd - ccad + ddac);
  MPoly<3> t1_rs_q = 4.0 * (+ bbcd - ccbd + ddbc);

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

  MPoly<3> eefg = e * e * f * g;
  MPoly<3> eefh = e * e * f * h;
  MPoly<3> eegh = e * e * g * h;
  MPoly<3> ffeg = f * f * e * g;
  MPoly<3> ffeh = f * f * e * h;
  MPoly<3> ffgh = f * f * g * h;
  MPoly<3> ggef = g * g * e * f;
  MPoly<3> ggeh = g * g * e * h;
  MPoly<3> ggfh = g * g * f * h;
  MPoly<3> hhef = h * h * e * f;
  MPoly<3> hheg = h * h * e * g;
  MPoly<3> hhfg = h * h * f * g;

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

#if 0
  printf("\tqr->p\tqr->s\tqs->p\tqs->r\trs->p\trs->q\n");
  printf("\tpq->r\tpq->s\tpr->q\tpr->s\tps->q\tps->r\n");
  printf("aabc\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(aabc, qr), p - qr),
         dot_product(gradient(aabc, qr), s - qr),
         dot_product(gradient(aabc, qs), p - qs),
         dot_product(gradient(aabc, qs), r - qs),
         dot_product(gradient(aabc, rs), p - rs),
         dot_product(gradient(aabc, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(aabc, pq), r - pq),
         dot_product(gradient(aabc, pq), s - pq),
         dot_product(gradient(aabc, pr), q - pr),
         dot_product(gradient(aabc, pr), s - pr),
         dot_product(gradient(aabc, ps), q - ps),
         dot_product(gradient(aabc, ps), r - ps));
  printf("aabd\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(aabd, qr), p - qr),
         dot_product(gradient(aabd, qr), s - qr),
         dot_product(gradient(aabd, qs), p - qs),
         dot_product(gradient(aabd, qs), r - qs),
         dot_product(gradient(aabd, rs), p - rs),
         dot_product(gradient(aabd, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(aabd, pq), r - pq),
         dot_product(gradient(aabd, pq), s - pq),
         dot_product(gradient(aabd, pr), q - pr),
         dot_product(gradient(aabd, pr), s - pr),
         dot_product(gradient(aabd, ps), q - ps),
         dot_product(gradient(aabd, ps), r - ps));
  printf("aacd\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(aacd, qr), p - qr),
         dot_product(gradient(aacd, qr), s - qr),
         dot_product(gradient(aacd, qs), p - qs),
         dot_product(gradient(aacd, qs), r - qs),
         dot_product(gradient(aacd, rs), p - rs),
         dot_product(gradient(aacd, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(aacd, pq), r - pq),
         dot_product(gradient(aacd, pq), s - pq),
         dot_product(gradient(aacd, pr), q - pr),
         dot_product(gradient(aacd, pr), s - pr),
         dot_product(gradient(aacd, ps), q - ps),
         dot_product(gradient(aacd, ps), r - ps));
  printf("bbac\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(bbac, qr), p - qr),
         dot_product(gradient(bbac, qr), s - qr),
         dot_product(gradient(bbac, qs), p - qs),
         dot_product(gradient(bbac, qs), r - qs),
         dot_product(gradient(bbac, rs), p - rs),
         dot_product(gradient(bbac, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(bbac, pq), r - pq),
         dot_product(gradient(bbac, pq), s - pq),
         dot_product(gradient(bbac, pr), q - pr),
         dot_product(gradient(bbac, pr), s - pr),
         dot_product(gradient(bbac, ps), q - ps),
         dot_product(gradient(bbac, ps), r - ps));
  printf("bbad\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(bbad, qr), p - qr),
         dot_product(gradient(bbad, qr), s - qr),
         dot_product(gradient(bbad, qs), p - qs),
         dot_product(gradient(bbad, qs), r - qs),
         dot_product(gradient(bbad, rs), p - rs),
         dot_product(gradient(bbad, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(bbad, pq), r - pq),
         dot_product(gradient(bbad, pq), s - pq),
         dot_product(gradient(bbad, pr), q - pr),
         dot_product(gradient(bbad, pr), s - pr),
         dot_product(gradient(bbad, ps), q - ps),
         dot_product(gradient(bbad, ps), r - ps));
  printf("bbcd\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(bbcd, qr), p - qr),
         dot_product(gradient(bbcd, qr), s - qr),
         dot_product(gradient(bbcd, qs), p - qs),
         dot_product(gradient(bbcd, qs), r - qs),
         dot_product(gradient(bbcd, rs), p - rs),
         dot_product(gradient(bbcd, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(bbcd, pq), r - pq),
         dot_product(gradient(bbcd, pq), s - pq),
         dot_product(gradient(bbcd, pr), q - pr),
         dot_product(gradient(bbcd, pr), s - pr),
         dot_product(gradient(bbcd, ps), q - ps),
         dot_product(gradient(bbcd, ps), r - ps));
  printf("ccab\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ccab, qr), p - qr),
         dot_product(gradient(ccab, qr), s - qr),
         dot_product(gradient(ccab, qs), p - qs),
         dot_product(gradient(ccab, qs), r - qs),
         dot_product(gradient(ccab, rs), p - rs),
         dot_product(gradient(ccab, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ccab, pq), r - pq),
         dot_product(gradient(ccab, pq), s - pq),
         dot_product(gradient(ccab, pr), q - pr),
         dot_product(gradient(ccab, pr), s - pr),
         dot_product(gradient(ccab, ps), q - ps),
         dot_product(gradient(ccab, ps), r - ps));
  printf("ccad\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ccad, qr), p - qr),
         dot_product(gradient(ccad, qr), s - qr),
         dot_product(gradient(ccad, qs), p - qs),
         dot_product(gradient(ccad, qs), r - qs),
         dot_product(gradient(ccad, rs), p - rs),
         dot_product(gradient(ccad, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ccad, pq), r - pq),
         dot_product(gradient(ccad, pq), s - pq),
         dot_product(gradient(ccad, pr), q - pr),
         dot_product(gradient(ccad, pr), s - pr),
         dot_product(gradient(ccad, ps), q - ps),
         dot_product(gradient(ccad, ps), r - ps));
  printf("ccbd\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ccbd, qr), p - qr),
         dot_product(gradient(ccbd, qr), s - qr),
         dot_product(gradient(ccbd, qs), p - qs),
         dot_product(gradient(ccbd, qs), r - qs),
         dot_product(gradient(ccbd, rs), p - rs),
         dot_product(gradient(ccbd, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ccbd, pq), r - pq),
         dot_product(gradient(ccbd, pq), s - pq),
         dot_product(gradient(ccbd, pr), q - pr),
         dot_product(gradient(ccbd, pr), s - pr),
         dot_product(gradient(ccbd, ps), q - ps),
         dot_product(gradient(ccbd, ps), r - ps));
  printf("ddab\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ddab, qr), p - qr),
         dot_product(gradient(ddab, qr), s - qr),
         dot_product(gradient(ddab, qs), p - qs),
         dot_product(gradient(ddab, qs), r - qs),
         dot_product(gradient(ddab, rs), p - rs),
         dot_product(gradient(ddab, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ddab, pq), r - pq),
         dot_product(gradient(ddab, pq), s - pq),
         dot_product(gradient(ddab, pr), q - pr),
         dot_product(gradient(ddab, pr), s - pr),
         dot_product(gradient(ddab, ps), q - ps),
         dot_product(gradient(ddab, ps), r - ps));
  printf("ddac\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ddac, qr), p - qr),
         dot_product(gradient(ddac, qr), s - qr),
         dot_product(gradient(ddac, qs), p - qs),
         dot_product(gradient(ddac, qs), r - qs),
         dot_product(gradient(ddac, rs), p - rs),
         dot_product(gradient(ddac, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ddac, pq), r - pq),
         dot_product(gradient(ddac, pq), s - pq),
         dot_product(gradient(ddac, pr), q - pr),
         dot_product(gradient(ddac, pr), s - pr),
         dot_product(gradient(ddac, ps), q - ps),
         dot_product(gradient(ddac, ps), r - ps));
  printf("ddbc\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ddbc, qr), p - qr),
         dot_product(gradient(ddbc, qr), s - qr),
         dot_product(gradient(ddbc, qs), p - qs),
         dot_product(gradient(ddbc, qs), r - qs),
         dot_product(gradient(ddbc, rs), p - rs),
         dot_product(gradient(ddbc, rs), q - rs));
  printf("\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\t%+.4f\n",
         dot_product(gradient(ddbc, pq), r - pq),
         dot_product(gradient(ddbc, pq), s - pq),
         dot_product(gradient(ddbc, pr), q - pr),
         dot_product(gradient(ddbc, pr), s - pr),
         dot_product(gradient(ddbc, ps), q - ps),
         dot_product(gradient(ddbc, ps), r - ps));
#endif

  should(vector_equal(gradient(quartic1, qr), gradient(quartic2, qr)));
  should(vector_equal(gradient(quartic1, qs), gradient(quartic2, qs)));
  should(vector_equal(gradient(quartic1, rs), gradient(quartic2, rs)));
}
