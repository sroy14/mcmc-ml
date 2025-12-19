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

#ifndef DATA_H
#define DATA_H

#include <o2scl/table.h>

#include "settings.h"

namespace mc2ml {

  class data {

  public:

    static const int ix_success=0;
    static const int ix_wgt_zero=1;
    static const int ix_grad_failure=2;

    o2scl::uniform_grid<double> m_grid;
    o2scl::table<> wgt_star;
    o2scl::table<> wgt_grid;
    
    std::vector<std::string> s_names;
    std::vector<std::string> s_units;
    std::vector<double> s_mass;
    std::vector<double> c_68;
    std::vector<double> d_68;

    // Number of stars in LMXB
    static const size_t n_lmxb=5;

    // Number of stars in HMXB
    static const size_t n_hmxb=11;

    // Number of stars in NS-NS
    static const size_t n_nsns=22;

    // Number of stars in NS-WD
    static const size_t n_nswd=32;

    // Total number of stars
    size_t n_stars;

    data() {};

    virtual ~data()=default;

    data(const data &dat) {
      m_grid=dat.m_grid;
      s_names=dat.s_names;
      s_units=dat.s_units;
      s_mass=dat.s_mass;
      c_68=dat.c_68;
      d_68=dat.d_68;
      wgt_star=dat.wgt_star;
      wgt_grid=dat.wgt_grid;
      n_stars=dat.n_stars;
    }

    data &operator=(const data &dat) {
      if (this!=&dat) {
        m_grid=dat.m_grid;
        s_names=dat.s_names;
        s_units=dat.s_units;
        s_mass=dat.s_mass;
        c_68=dat.c_68;
        d_68=dat.d_68;
        wgt_star=dat.wgt_star;
        wgt_grid=dat.wgt_grid;
        n_stars=dat.n_stars;
      }
      return *this;
    }

    virtual void get_param_info(std::vector<std::string> &,
                                std::vector<std::string> &,
                                std::vector<double> &,
                                std::vector<double> &,
                                std::shared_ptr<settings>);

    virtual void set_init_point(std::vector<double> &,
                                std::shared_ptr<settings>);

    virtual void load_data(std::shared_ptr<settings>);

  };


  class solver {

  public:

    double funct_1d(double, double &, double &);
    double solve_1d(double, double);

  };

}

#endif // DATA_H