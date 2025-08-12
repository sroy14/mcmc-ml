/* ==============================================================
   pdf_demo.cpp  —  compile with:  g++ -std=c++17 -O3 pdf_demo.cpp -o demo
   ============================================================== */
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <o2scl/constants.h>   // provides o2scl_const::pi

// ------------------------------------------------------------------
//  Renamed PDF class (exactly your latest version)
// ------------------------------------------------------------------
class pdf {

public:

    static inline double s_norm(double u) {
        return std::exp(-0.5*u*u) /
               std::sqrt(2.0 * o2scl_const::pi);
    }

    static inline double c_norm(double u) {
        return 0.5 * (1.0 + std::erf(u / std::sqrt(2.0)));
    }

    // ------ skew‑normal -------------------------------------------------
    static inline double skewed_norm(double x,double m,double s,double a) {
        double z = (x - m) / s;
        return 2.0 / s * s_norm(z) * c_norm(a * z);
    }

    static inline double dsn_dx(double x,double m,double s,double a) {
        double z = (x - m) / s;
        return 2.0 / (s * s) * s_norm(z)
             * ( -z * c_norm(a * z) + a * s_norm(a * z) );
    }

    static inline double dsn_dm(double x,double m,double s,double a) {
        return -dsn_dx(x, m, s, a);
    }

    static inline double dsn_ds(double x,double m,double s,double a) {
        double z = (x - m) / s;
        return 2.0 / (s * s) * s_norm(z)
             * ( (z * z - 1.0) * c_norm(a * z) - a * z * s_norm(a * z) );
    }

    static inline double dsn_da(double x,double m,double s,double a) {
        double z = (x - m) / s;
        return 2.0 / (s * s) * (x - m) * s_norm(z) * s_norm(a * z);
    }

    // ------ asymmetric‑normal ------------------------------------------
    static inline double asym_norm(double x,double c,double d) {
        double k = 2.0 / (d * (c + 1.0 / c));
        double u = (x < 0.0) ? (c * x / d) : (x / (c * d));
        return k * s_norm(u);
    }

    static inline double dan_dx(double x,double c,double d) {
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
//  Five‑point central difference for high accuracy
// ------------------------------------------------------------------
double cdiff5(const std::function<double(double)>& f, double x)
{
    const double eps  = 2.2204460492503131e-16;      //  double ε
    double h = std::pow(eps, 0.2) *
               (std::fabs(x) > 1.0 ? std::fabs(x) : 1.0);

    return ( -f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h) ) / (12*h);
}

// ------------------------------------------------------------------
//  Pretty printer
// ------------------------------------------------------------------
void print_cmp(const char* tag, double ana, double num)
{
    std::cout << "  " << tag
              << "  analytic " << std::setprecision(16) << std::scientific << ana
              << "  numeric "  << num
              << "  abs‑diff " << std::fabs(ana - num) << '\n';
}

// ------------------------------------------------------------------
//  main() – run one test point
// ------------------------------------------------------------------
int main()
{
    double x = 0.7,  m = 0.2,  s = 1.3,  a = 0.8;
    double c = 1.5,  d = 0.6;

    std::cout << "\nSkewed‑normal\n";
    auto f_sn = [&](double xx){ return pdf::skewed_norm(xx, m, s, a); };

    print_cmp("dx", pdf::dsn_dx(x,m,s,a), cdiff5(f_sn, x));
    print_cmp("dm", pdf::dsn_dm(x,m,s,a),
              cdiff5([&](double mm){ return pdf::skewed_norm(x,mm,s,a); }, m));
    print_cmp("ds", pdf::dsn_ds(x,m,s,a),
              cdiff5([&](double ss){ return pdf::skewed_norm(x,m,ss,a); }, s));
    print_cmp("da", pdf::dsn_da(x,m,s,a),
              cdiff5([&](double aa){ return pdf::skewed_norm(x,m,s,aa); }, a));

    std::cout << "\nAsymmetric‑normal  x = -0.8\n";
    auto f_an = [&](double xx){ return pdf::asym_norm(xx,c,d); };
    print_cmp("dx", pdf::dan_dx(-0.8,c,d), cdiff5(f_an, -0.8));

    std::cout << "\nAsymmetric‑normal  x =  0.8\n";
    print_cmp("dx", pdf::dan_dx( 0.8,c,d), cdiff5(f_an,  0.8));
}
