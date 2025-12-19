// test.cpp  — gradient of product likelihood using uBLAS vectors
// Build demo: g++ -std=c++17 -O3 -DTEST_MAIN test.cpp -o test && ./test
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include <o2scl/constants.h>                         // o2scl_const::pi
#include <boost/numeric/ublas/vector.hpp>

using ubvector = boost::numeric::ublas::vector<double>;

// New callback signature (as provided)
struct data; // fwd
typedef std::function<int(std::size_t, const ubvector&, double&, data&)> point_funct;

// ------------------------------------------------------------------
// Minimal data holder (sizes inferred at runtime)
//   s_mass[i][j] : measured mass m_{ij}
//   c_68[i][j]   : asym-normal c_{ij}
//   d_68[i][j]   : asym-normal d_{ij}
// ------------------------------------------------------------------
struct data {
  std::vector<std::vector<double>> s_mass;
  std::vector<std::vector<double>> c_68;
  std::vector<std::vector<double>> d_68;
};

// ------------------------------------------------------------------
// pdf utilities (validated formulas)
// ------------------------------------------------------------------
struct pdf {

  static inline double s_norm(double u) {
    return std::exp(-0.5*u*u) / std::sqrt(2.0 * o2scl_const::pi);
  }
  static inline double c_norm(double u) {
    return 0.5 * (1.0 + std::erf(u / std::sqrt(2.0)));
  }

  // Skewed normal: SN(x; m,s,a) = 2/s * φ((x-m)/s) * Φ(a*(x-m)/s)
  static inline double skewed_norm(double x, double m, double s, double a) {
    double z = (x - m) / s;
    return 2.0 / s * s_norm(z) * c_norm(a * z);
  }

  // Asymmetric normal: AN(w; c,d) with piecewise stretch
  static inline double asym_norm(double x, double c, double d) {
    double k = 2.0 / (d * (c + 1.0 / c));
    double u = (x < 0.0) ? (c * x / d) : (x / (c * d));
    return k * s_norm(u);
  }

  // ∂SN/∂x (needed for M-derivative)
  static inline double dsn_dx(double x, double m, double s, double a) {
    double z = (x - m) / s;
    return 2.0 / (s * s) * s_norm(z) *
           (-z * c_norm(a * z) + a * s_norm(a * z));
  }

  // ∂AN/∂x (x is w); used via chain rule for M (w = m − M ⇒ ∂/∂M = −∂/∂w)
  static inline double dan_dx(double x, double c, double d) {
    double k = 2.0 / (d * (c + 1.0 / c));
    if (x < 0.0) {
      double u = c * x / d;
      return -k * c * c * x / (d * d) * s_norm(u);
    } else {
      double u = x / (c * d);
      return -k * x / (c * c * d * d) * s_norm(u);
    }
  }
};

// ------------------------------------------------------------------
// Likelihood L(pars; d) with parameter packing:
//   pars = [ mu[0..P-1], sig[0..P-1], alp[0..P-1], M_ij (i-major, j-minor) ]
// Sizes:  P = d.s_mass.size(),  N_i = d.s_mass[i].size()
// ------------------------------------------------------------------
static double likelihood(const ubvector &pars, data &d)
{
  const int P = static_cast<int>(d.s_mass.size());
  std::vector<int> Ni(P);
  long long sumN = 0;
  for (int i=0;i<P;++i){ Ni[i] = static_cast<int>(d.s_mass[i].size()); sumN += Ni[i]; }

  const long long off_mu  = 0;
  const long long off_sig = off_mu  + P;
  const long long off_alp = off_sig + P;
  const long long off_M   = off_alp + P;

  double L = 1.0;
  long long base = off_M;
  for (int i=0;i<P;++i){
    for (int j=0;j<Ni[i];++j){
      const double mu  = pars(off_mu  + i);
      const double sg  = pars(off_sig + i);
      const double al  = pars(off_alp + i);
      const double Mij = pars(base + j);

      const double SN = pdf::skewed_norm(Mij, mu, sg, al);

      const double w   = d.s_mass[i][j] - Mij;   // w = m − M
      const double cij = d.c_68[i][j];
      const double dij = d.d_68[i][j];
      const double AN  = pdf::asym_norm(w, cij, dij);

      L *= AN * SN;
    }
    base += Ni[i];
  }
  return L;
}

// ------------------------------------------------------------------
// 1) Analytic gradient: grad matches the parameter packing of 'pars'.
//    Signature per request: int deriv(const ubvector&, point_funct&, ubvector&, data&)
//    pf is ignored (present for interface compatibility).
// ------------------------------------------------------------------
int deriv(const ubvector &pars, point_funct &pf, ubvector &grad, data &d) {

  size_t P = 1;
  std::vector<int> Ni(P);
  long long sumN = 0;
  for (int i=0;i<P;++i){ 
    Ni[i] = static_cast<int>(d.s_mass[i].size()); 
    sumN += Ni[i]; 
  }

  const long long n_expected = 3LL*P + sumN;
  if (static_cast<long long>(pars.size()) != n_expected) {
    std::cerr << "deriv(): parameter size mismatch. Got " << pars.size()
              << ", expected " << n_expected << "\n";
    return 1;
  }

  grad.resize(pars.size(), false);
  for (std::size_t k=0;k<grad.size();++k) grad(k)=0.0;

  const long long off_mu  = 0;
  const long long off_sig = off_mu  + P;
  const long long off_alp = off_sig + P;
  const long long off_M   = off_alp + P;

  // Precompute likelihood once (used to scale log-derivatives)
  const double L = likelihood(pars, d);

  // Prefix for indexing the M-block
  std::vector<long long> pref(P+1,0);
  for (int i=0;i<P;++i) pref[i+1] = pref[i] + Ni[i];

  // Accumulate population-wise sums for μ,σ,α
  std::vector<double> sum_mu(P,0.0), sum_sig(P,0.0), sum_alp(P,0.0);

  // Loop over all (i,j)
  for (int i=0;i<P;++i){
    const double mu  = pars(off_mu  + i);
    const double sg  = pars(off_sig + i);
    const double al  = pars(off_alp + i);

    for (int j=0;j<Ni[i];++j){
      const long long idxM = off_M + pref[i] + j;
      const double Mij = pars(idxM);

      // SN pieces
      const double z   = (Mij - mu) / sg;
      const double Phi = pdf::c_norm(al * z);
      const double phi = pdf::s_norm(al * z);
      const double ratio = phi / Phi;

      // log-derivatives from SN
      const double dlogSN_dmu  = (1.0/sg) * (  z - al * ratio );
      const double dlogSN_dsig = (1.0/sg) * ( (z*z - 1.0) - al * z * ratio );
      const double dlogSN_dalp =              z * ratio;
      const double dlogSN_dM   = (1.0/sg) * ( -z + al * ratio );

      sum_mu[i]  += dlogSN_dmu;
      sum_sig[i] += dlogSN_dsig;
      sum_alp[i] += dlogSN_dalp;

      // AN contribution via w = m − M
      const double w   = d.s_mass[i][j] - Mij;
      const double cij = d.c_68[i][j];
      const double dij = d.d_68[i][j];
      const double AN  = pdf::asym_norm(w, cij, dij);
      const double dANdM_over_AN = (-pdf::dan_dx(w, cij, dij)) / AN; // chain rule

      // Full derivative for M_ij: L * ( dlogSN/dM + dlogAN/dM )
      grad(idxM) = L * ( dlogSN_dM + dANdM_over_AN );
    }
  }

  // Scale population sums by L to get ∂L/∂μ_i, ∂L/∂σ_i, ∂L/∂α_i
  for (int i=0;i<P;++i){
    grad(off_mu  + i) = L * sum_mu[i];
    grad(off_sig + i) = L * sum_sig[i];
    grad(off_alp + i) = L * sum_alp[i];
  }

  return 0;
}

// ------------------------------------------------------------------
// 2) Numeric gradient via 5-point stencil, component-wise.
//    Signature per request: int num_deriv(const ubvector&, point_funct&, ubvector&, data&)
//    Uses 'f' to evaluate the scalar function value; returns 0 on success.
// ------------------------------------------------------------------
int num_deriv(const ubvector &x, point_funct &f, ubvector &g, data &dat)
{
  const std::size_t n = x.size();
  g.resize(n, false);

  const double eps = 2.2204460492503131e-16; // double ε
  for (std::size_t k=0;k<n;++k){
    const double h = std::pow(eps, 0.2) * (std::fabs(x(k)) > 1.0 ? std::fabs(x(k)) : 1.0);

    ubvector xp2 = x, xp = x, xm = x, xm2 = x;
    xp2(k) = x(k) + 2.0*h;
    xp (k) = x(k) + h;
    xm (k) = x(k) - h;
    xm2(k) = x(k) - 2.0*h;

    double fp2=0.0, fp=0.0, fm=0.0, fm2=0.0;

    // size_t argument is unused by our scalar objective; pass 0
    int rc = 0;
    rc |= f(0, xp2, fp2, dat);
    rc |= f(0, xp , fp , dat);
    rc |= f(0, xm , fm , dat);
    rc |= f(0, xm2, fm2, dat);
    if (rc) return rc;

    g(k) = ( -fp2 + 8.0*fp - 8.0*fm + fm2 ) / (12.0*h);
  }
  return 0;
}

// ------------------------------------------------------------------
// Optional demonstration (enable with -DTEST_MAIN)
// ------------------------------------------------------------------
#ifdef TEST_MAIN
// Wrap our internal likelihood into the required point_funct signature
static int like_pf(std::size_t /*unused*/, const ubvector &p, double &val, data &d)
{
  val = likelihood(p, d);
  return 0;
}

int main(){
  // Tiny P=1 example (user normally supplies these externally)
  data d;
  d.s_mass = {{0.7, -0.3}};   // m_{1j}
  d.c_68   = {{1.5,  1.2}};
  d.d_68   = {{0.6,  0.7}};

  // pars = [mu0, sig0, alp0, M_10, M_11]
  ubvector pars(5);
  pars(0)=0.1;  pars(1)=1.2;  pars(2)=0.7;  pars(3)=0.7;  pars(4)=-0.3;

  point_funct pf = like_pf;

  ubvector ga, gn;
  int rc1 = deriv(pars, pf, ga, d);
  int rc2 = num_deriv(pars, pf, gn, d);
  if (rc1 || rc2) { std::cerr<<"error: rc="<<(rc1?rc1:rc2)<<"\n"; return 1; }

  std::cout.setf(std::ios::scientific); std::cout<<std::setprecision(16);
  for (std::size_t k=0;k<pars.size();++k){
    std::cout<<"k="<<k<<"  ana="<<ga(k)<<"  num="<<gn(k)
             <<"  abs-diff="<<std::fabs(ga(k)-gn(k))<<"\n";
  }
  return 0;
}
#endif
