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

    o2scl::table<> grid;
    o2scl::uniform_grid<double> m_grid;

    std::vector<double> m_dt;
    std::vector<double> c_68;
    std::vector<double> d_68;
    static const size_t n_stars=5;

    data() {
    }

    data(const data &dat) {
      grid=dat.grid;
    }

    data &operator=(const data &dat) {
      if (this!=&dat) {
        grid=dat.grid;
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

    virtual void load_mass_data();

  };


  class solver {

  public:

    double funct_x(double, double &, double &);
    double solve_x(double, double);

  };

}

#endif // DATA_H