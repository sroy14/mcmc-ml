/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2025, Mahmudul Hasan Anik, Satyajit Roy, and 
  Andrew W. Steiner
  
  This file is part of mcmc-ml.

  Mcmc-ml is free software; you can redistribute it and/or modify it 
  under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  Mcmc-ml is distributed in the hope that it will be useful, but 
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with mcmc-ml. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/

#ifndef BASE_H
#define BASE_H

#include <o2scl/constants.h>

#include "data.h"
#include "settings.h"

namespace mc2ml {

  typedef boost::numeric::ublas::vector<double> ubvector;
  
  typedef std::function<int(size_t, const ubvector &, 
                            double &, data &)> point_funct;
  
  typedef std::function<int(const ubvector &, double, 
                            std::vector<double> &, data &)> fill_funct;

  typedef std::function<int(size_t, const ubvector &, point_funct &, 
                            ubvector &, data &, bool &)> deriv_funct;
  
  class base {
  
  public:
    
    int n_threads;
    
    std::shared_ptr<data> dat;
    std::shared_ptr<settings> set;
    
    std::string mc_type;
    std::string ml_type;

    o2scl::vec_index pvi;

    base() {
      n_threads=1;
    }

    virtual ~base() {
    }

    virtual int point(const ubvector &, std::ofstream &, double &, data &);
    virtual int deriv(const ubvector &, point_funct &, ubvector &, data &, bool &);
    virtual int fill(const ubvector &, double, std::vector<double> &, data &);

  };


  class pdf {

  public:

    static inline double s_norm(double u) {
      return std::exp(-0.5*u*u)/std::sqrt(2.0*o2scl_const::pi);
    }

    static inline double c_norm(double u) {
      return 0.5*(1.0+std::erf(u/std::sqrt(2.0)));
    }

    static inline double skewed_norm(double x, double m, double s, double a) {
      double z=(x-m)/s;
      return 2.0/s*s_norm(z)*c_norm(a*z);
    }

    static inline double dsn_dx(double x, double m, double s, double a) {
      double z=(x-m)/s;
      return 2.0/(s*s)*s_norm(z)*(-z*c_norm(a*z)+a*s_norm(a*z));
    }

    static inline double dsn_dm(double x, double m, double s, double a) {
      return -dsn_dx(x, m, s, a);
    }

    static inline double dsn_ds(double x, double m, double s, double a) {
      double z=(x-m)/s;
      return 2.0/(s*s)*s_norm(z)*((z*z-1.0)*c_norm(a*z)-a*z*s_norm(a*z));
    }

    static inline double dsn_da(double x, double m, double s, double a) {
      double z=(x-m)/s;
      return 2.0/(s*s)*(x-m)*s_norm(z)*s_norm(a*z);
    }

    static inline double asym_norm(double x, double c, double d) {
      double k=2.0/(d*(c+1.0/c));
      double u;
      if (x<0.0) u=c*x/d;
      else u=x/(c*d);
      return k*s_norm(u);
    }

    static inline double dan_dx(double x, double c, double d) {
      double k=2.0/(d*(c+1.0/c));
      if (x<0.0) {
          double u=c*x/d;
          return -k*c*c*x/(d*d)*s_norm(u);
      } else {
          double u=x/(c*d);
          return -k*x/(c*c*d*d)*s_norm(u);
      }
    }
  };



}

#endif // MAIN_H