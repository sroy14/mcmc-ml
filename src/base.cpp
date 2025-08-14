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

#include "base.h"
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace mc2ml;


int base::fill(const ubvector &pars, double wgt, 
               vector<double> &line, data &d) {
  
  if (set->inc_lmxb) {
    for (size_t i=0; i<d.n_stars; i++) {
      line.push_back(d.wgt_star.get("wgt", i));
    }
  } else {
    for (size_t i=0; i<set->grid_size; i++) {
      line.push_back(d.wgt_grid.get("wgt", i));
    }
  }

  return 0;

}


int base::point(const ubvector &pars, std::ofstream &sout, 
                double &log_wgt, data &d) {
  
  double mean=pars[pvi["mean"]];
  double width=pow(10.0, pars[pvi["log10_std"]]);
  double skew=pars[pvi["skewness"]];
  log_wgt=0.0;

  if (set->inc_lmxb) {

    if (d.wgt_star.get_ncolumns()==0) {
      d.wgt_star.new_column("wgt");
      d.wgt_star.set_nlines(d.n_stars);
    }

    for (size_t i=0; i<d.n_stars; i++) {

      double m_dat=d.s_mass[i];
      double asym=d.c_68[i];
      double scale=d.d_68[i];
      double m_par=pars[3+i];
      double wgt=pdf::asym_norm(m_dat-m_par, asym, scale) * 
                 pdf::skewed_norm(m_par, mean, width, skew);

      if (wgt<=0.0) {
        sout << "base::point(): LMXB star " << i 
             << " returned zero weight." << endl;
        // log_wgt=double(d.ix_wgt_zero)-100.0;
        return d.ix_wgt_zero;
      }

      d.wgt_star.set("wgt", i, wgt);
      log_wgt+=log(wgt);

    }

  } else { // !set->inc_lmxb

    if (d.wgt_grid.get_ncolumns()==0) {
      d.wgt_grid.new_column("wgt");
      d.wgt_grid.set_nlines(set->grid_size);
    }

    for (size_t i=0; i<set->grid_size; i++) {
      double m_val=d.m_grid[i];
      double wgt=pdf::skewed_norm(m_val, mean, width, skew);

      if (wgt<=0.0) {
        sout << "base::point(): grid point " << i 
             << " returned zero weight." << endl;
        // log_wgt=double(d.ix_wgt_zero)-100.0;
        return d.ix_wgt_zero;
      }

      d.wgt_grid.set("wgt", i, wgt);
      log_wgt+=log(wgt);
    }
  }

  return 0;

}

int base::deriv(const ubvector &pars, point_funct &pf, 
                ubvector &grad, data &d, bool &success) {
  return 0;
}