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

#include "main.h"
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace mc2ml;


int main::fill(const ubvector &pars, double wgt, 
               vector<double> &line, data &dat) {
  return 0;
}


int main::point(const ubvector &pars, std::ofstream &sout, 
                double &log_wgt, data &dat) {
  
  double mean=pars[pvi["mean"]];
  double width=pow(10.0, pars[pvi["log10_var"]]);
  double skew=pars[pvi["skewness"]];
  log_wgt=0.0;

  if (set->inc_lmxb) {
    for (size_t i=0; i<dat.n_stars; i++) {
      double m_dat=dat.m_dt[i];
      double asym=dat.c_68[i];
      double scale=dat.d_68[i];
      double m_par=pars[3+i];
      double wgt=pdf::asym_norm(m_dat-m_par, asym, scale) * 
                 pdf::skewed_norm(m_par, mean, width, skew);
      if (wgt<=0.0) {
        sout << "main::point(): LMXB star " << i 
             << " returned zero weight." << endl;
        // log_wgt=double(dat.ix_wgt_zero)-100.0;
        return dat.ix_wgt_zero;
      }
      log_wgt+=log(wgt);
    }
  } else {
    for (size_t i=0; i<set->grid_size; i++) {
      double m_val=dat.m_grid[i];
      double wgt=pdf::skewed_norm(m_val, mean, width, skew);
      if (wgt<=0.0) {
        sout << "main::point(): grid point " << i 
             << " returned zero weight." << endl;
        // log_wgt=double(dat.ix_wgt_zero)-100.0;
        return dat.ix_wgt_zero;
      }
      log_wgt+=log(wgt);
    }
  }

  return 0;

}

int main::deriv(const ubvector &pars, point_funct &pf, 
                ubvector &grad, data &dat, bool &success) {
  return 0;
}