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

  should(double_equal(aab(p), 0.0));
  should(double_equal(aab(q), 0.0));
  should(double_equal(aab(r), 0.0));
  should(double_equal(aab(s), 0.0));

  should(double_equal(aac(p), 0.0));
  should(double_equal(aac(q), 0.0));
  should(double_equal(aac(r), 0.0));
  should(double_equal(aac(s), 0.0));

  should(double_equal(aad(p), 0.0));
  should(double_equal(aad(q), 0.0));
  should(double_equal(aad(r), 0.0));
  should(double_equal(aad(s), 0.0));

  should(double_equal(bba(p), 0.0));
  should(double_equal(bba(q), 0.0));
  should(double_equal(bba(r), 0.0));
  should(double_equal(bba(s), 0.0));

  should(double_equal(bbc(p), 0.0));
  should(double_equal(bbc(q), 0.0));
  should(double_equal(bbc(r), 0.0));
  should(double_equal(bbc(s), 0.0));

  should(double_equal(bbd(p), 0.0));
  should(double_equal(bbd(q), 0.0));
  should(double_equal(bbd(r), 0.0));
  should(double_equal(bbd(s), 0.0));

  should(double_equal(cca(p), 0.0));
  should(double_equal(cca(q), 0.0));
  should(double_equal(cca(r), 0.0));
  should(double_equal(cca(s), 0.0));

  should(double_equal(ccb(p), 0.0));
  should(double_equal(ccb(q), 0.0));
  should(double_equal(ccb(r), 0.0));
  should(double_equal(ccb(s), 0.0));

  should(double_equal(ccd(p), 0.0));
  should(double_equal(ccd(q), 0.0));
  should(double_equal(ccd(r), 0.0));
  should(double_equal(ccd(s), 0.0));

  should(double_equal(dda(p), 0.0));
  should(double_equal(dda(q), 0.0));
  should(double_equal(dda(r), 0.0));
  should(double_equal(dda(s), 0.0));

  should(double_equal(ddb(p), 0.0));
  should(double_equal(ddb(q), 0.0));
  should(double_equal(ddb(r), 0.0));
  should(double_equal(ddb(s), 0.0));

  should(double_equal(ddc(p), 0.0));
  should(double_equal(ddc(q), 0.0));
  should(double_equal(ddc(r), 0.0));
  should(double_equal(ddc(s), 0.0));

  should(double_equal(dot_product(gradient(aab, p), q - p), 1.0));
  should(double_equal(dot_product(gradient(aab, p), r - p), 0.0));
  should(double_equal(dot_product(gradient(aab, p), s - p), 0.0));
  should(double_equal(dot_product(gradient(aac, p), q - p), 0.0));
  should(double_equal(dot_product(gradient(aac, p), r - p), 1.0));
  should(double_equal(dot_product(gradient(aac, p), s - p), 0.0));
  should(double_equal(dot_product(gradient(aad, p), q - p), 0.0));
  should(double_equal(dot_product(gradient(aad, p), r - p), 0.0));
  should(double_equal(dot_product(gradient(aad, p), s - p), 1.0));

  should(double_equal(dot_product(gradient(bba, q), p - q), 1.0));
  should(double_equal(dot_product(gradient(bba, q), r - q), 0.0));
  should(double_equal(dot_product(gradient(bba, q), s - q), 0.0));
  should(double_equal(dot_product(gradient(bbc, q), p - q), 0.0));
  should(double_equal(dot_product(gradient(bbc, q), r - q), 1.0));
  should(double_equal(dot_product(gradient(bbc, q), s - q), 0.0));
  should(double_equal(dot_product(gradient(bbd, q), p - q), 0.0));
  should(double_equal(dot_product(gradient(bbd, q), r - q), 0.0));
  should(double_equal(dot_product(gradient(bbd, q), s - q), 1.0));

  should(double_equal(dot_product(gradient(cca, r), p - r), 1.0));
  should(double_equal(dot_product(gradient(cca, r), q - r), 0.0));
  should(double_equal(dot_product(gradient(cca, r), s - r), 0.0));
  should(double_equal(dot_product(gradient(ccb, r), p - r), 0.0));
  should(double_equal(dot_product(gradient(ccb, r), q - r), 1.0));
  should(double_equal(dot_product(gradient(ccb, r), s - r), 0.0));
  should(double_equal(dot_product(gradient(ccd, r), p - r), 0.0));
  should(double_equal(dot_product(gradient(ccd, r), q - r), 0.0));
  should(double_equal(dot_product(gradient(ccd, r), s - r), 1.0));

  should(double_equal(dot_product(gradient(dda, s), p - s), 1.0));
  should(double_equal(dot_product(gradient(dda, s), q - s), 0.0));
  should(double_equal(dot_product(gradient(dda, s), r - s), 0.0));
  should(double_equal(dot_product(gradient(ddb, s), p - s), 0.0));
  should(double_equal(dot_product(gradient(ddb, s), q - s), 1.0));
  should(double_equal(dot_product(gradient(ddb, s), r - s), 0.0));
  should(double_equal(dot_product(gradient(ddc, s), p - s), 0.0));
  should(double_equal(dot_product(gradient(ddc, s), q - s), 0.0));
  should(double_equal(dot_product(gradient(ddc, s), r - s), 1.0));

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

  should(double_equal(eef(t), 0.0));
  should(double_equal(eef(q), 0.0));
  should(double_equal(eef(r), 0.0));
  should(double_equal(eef(s), 0.0));

  should(double_equal(eeg(t), 0.0));
  should(double_equal(eeg(q), 0.0));
  should(double_equal(eeg(r), 0.0));
  should(double_equal(eeg(s), 0.0));

  should(double_equal(eeh(t), 0.0));
  should(double_equal(eeh(q), 0.0));
  should(double_equal(eeh(r), 0.0));
  should(double_equal(eeh(s), 0.0));

  should(double_equal(ffe(t), 0.0));
  should(double_equal(ffe(q), 0.0));
  should(double_equal(ffe(r), 0.0));
  should(double_equal(ffe(s), 0.0));

  should(double_equal(ffg(t), 0.0));
  should(double_equal(ffg(q), 0.0));
  should(double_equal(ffg(r), 0.0));
  should(double_equal(ffg(s), 0.0));

  should(double_equal(ffh(t), 0.0));
  should(double_equal(ffh(q), 0.0));
  should(double_equal(ffh(r), 0.0));
  should(double_equal(ffh(s), 0.0));

  should(double_equal(gge(t), 0.0));
  should(double_equal(gge(q), 0.0));
  should(double_equal(gge(r), 0.0));
  should(double_equal(gge(s), 0.0));

  should(double_equal(ggf(t), 0.0));
  should(double_equal(ggf(q), 0.0));
  should(double_equal(ggf(r), 0.0));
  should(double_equal(ggf(s), 0.0));

  should(double_equal(ggh(t), 0.0));
  should(double_equal(ggh(q), 0.0));
  should(double_equal(ggh(r), 0.0));
  should(double_equal(ggh(s), 0.0));

  should(double_equal(hhe(t), 0.0));
  should(double_equal(hhe(q), 0.0));
  should(double_equal(hhe(r), 0.0));
  should(double_equal(hhe(s), 0.0));

  should(double_equal(hhf(t), 0.0));
  should(double_equal(hhf(q), 0.0));
  should(double_equal(hhf(r), 0.0));
  should(double_equal(hhf(s), 0.0));

  should(double_equal(hhg(t), 0.0));
  should(double_equal(hhg(q), 0.0));
  should(double_equal(hhg(r), 0.0));
  should(double_equal(hhg(s), 0.0));

  should(double_equal(dot_product(gradient(eef, t), q - t), 1.0));
  should(double_equal(dot_product(gradient(eef, t), r - t), 0.0));
  should(double_equal(dot_product(gradient(eef, t), s - t), 0.0));
  should(double_equal(dot_product(gradient(eeg, t), q - t), 0.0));
  should(double_equal(dot_product(gradient(eeg, t), r - t), 1.0));
  should(double_equal(dot_product(gradient(eeg, t), s - t), 0.0));
  should(double_equal(dot_product(gradient(eeh, t), q - t), 0.0));
  should(double_equal(dot_product(gradient(eeh, t), r - t), 0.0));
  should(double_equal(dot_product(gradient(eeh, t), s - t), 1.0));

  should(double_equal(dot_product(gradient(ffe, q), t - q), 1.0));
  should(double_equal(dot_product(gradient(ffe, q), r - q), 0.0));
  should(double_equal(dot_product(gradient(ffe, q), s - q), 0.0));
  should(double_equal(dot_product(gradient(ffg, q), t - q), 0.0));
  should(double_equal(dot_product(gradient(ffg, q), r - q), 1.0));
  should(double_equal(dot_product(gradient(ffg, q), s - q), 0.0));
  should(double_equal(dot_product(gradient(ffh, q), t - q), 0.0));
  should(double_equal(dot_product(gradient(ffh, q), r - q), 0.0));
  should(double_equal(dot_product(gradient(ffh, q), s - q), 1.0));

  should(double_equal(dot_product(gradient(gge, r), t - r), 1.0));
  should(double_equal(dot_product(gradient(gge, r), q - r), 0.0));
  should(double_equal(dot_product(gradient(gge, r), s - r), 0.0));
  should(double_equal(dot_product(gradient(ggf, r), t - r), 0.0));
  should(double_equal(dot_product(gradient(ggf, r), q - r), 1.0));
  should(double_equal(dot_product(gradient(ggf, r), s - r), 0.0));
  should(double_equal(dot_product(gradient(ggh, r), t - r), 0.0));
  should(double_equal(dot_product(gradient(ggh, r), q - r), 0.0));
  should(double_equal(dot_product(gradient(ggh, r), s - r), 1.0));

  should(double_equal(dot_product(gradient(hhe, s), t - s), 1.0));
  should(double_equal(dot_product(gradient(hhe, s), q - s), 0.0));
  should(double_equal(dot_product(gradient(hhe, s), r - s), 0.0));
  should(double_equal(dot_product(gradient(hhf, s), t - s), 0.0));
  should(double_equal(dot_product(gradient(hhf, s), q - s), 1.0));
  should(double_equal(dot_product(gradient(hhf, s), r - s), 0.0));
  should(double_equal(dot_product(gradient(hhg, s), t - s), 0.0));
  should(double_equal(dot_product(gradient(hhg, s), q - s), 0.0));
  should(double_equal(dot_product(gradient(hhg, s), r - s), 1.0));
}
