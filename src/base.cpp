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
               vector<double> &line, data &dat) {
  
  if (set->inc_lmxb) {
    for (size_t i=0; i<dat.n_stars; i++) {
      line.push_back(dat.wgt_star.get("wgt", i));
    }
  } else {
    for (size_t i=0; i<set->grid_size; i++) {
      line.push_back(dat.wgt_grid.get("wgt", i));
    }
  }

  return 0;

}


int base::point(const ubvector &pars, std::ofstream &sout, 
                double &log_wgt, data &dat) {
  
  double mean=pars[pvi["mean"]];
  double width=pow(10.0, pars[pvi["log10_std"]]);
  double skew=pars[pvi["skewness"]];
  log_wgt=0.0;

  if (set->inc_lmxb) {

    if (dat.wgt_star.get_ncolumns()==0) {
      dat.wgt_star.set_nlines(dat.n_stars);
      dat.wgt_star.new_column("wgt");
    }

    for (size_t i=0; i<dat.n_stars; i++) {

      double m_dat=dat.s_mass[i];
      double asym=dat.c_68[i];
      double scale=dat.d_68[i];
      double m_par=pars[3+i];
      double wgt=pdf::asym_norm(m_dat-m_par, asym, scale) * 
                 pdf::skewed_norm(m_par, mean, width, skew);

      if (wgt<=0.0) {
        sout << "base::point(): LMXB star " << i 
             << " returned zero weight." << endl;
        // log_wgt=double(dat.ix_wgt_zero)-100.0;
        return dat.ix_wgt_zero;
      }

      dat.wgt_star.set("wgt", i, wgt);
      log_wgt+=log(wgt);

    }

  } else { // !set->inc_lmxb

    if (dat.wgt_grid.get_ncolumns()==0) {
      dat.wgt_grid.new_column("wgt");
      dat.wgt_grid.set_nlines(set->grid_size);
    }

    double sum=0.0, dx=dat.m_grid[1]-dat.m_grid[0];

    for (size_t i=0; i<set->grid_size; i++) {

      double m_val=dat.m_grid[i];
      double wgt=pdf::skewed_norm(m_val, mean, width, skew);

      if (wgt<=0.0) {
        sout << "base::point(): grid point " << i 
             << " returned zero weight." << endl;
        // log_wgt=double(dat.ix_wgt_zero)-100.0;
        return dat.ix_wgt_zero;
      }

      dat.wgt_grid.set("wgt", i, wgt);
      sum+=wgt;
      
    }

    double integral=sum*dx;
    log_wgt+=log(integral);
  
  }

  return 0;

}

int base::deriv(const ubvector &pars, point_funct &pf, 
                ubvector &grad, data &dat) {
  
  size_t n_params=pars.size();
  if (grad.size()!=n_params) grad.resize(n_params);

  double weights=0.0;

  for (size_t i=0; i<dat.n_stars; i++) {
    weights+=dat.wgt_star.get("wgt", i);
  }

  // Main loop over all (i,j)
  for (size_t i=0; i<dat.n_pops; i++) {
    
    double mean=pars[pvi["mean"]];
    double width=pow(10.0, pars[pvi["log10_std"]]);
    double skew=pars[pvi["skewness"]];

    for (size_t j=0; j<dat.n_lmxb; j++) {
      double m_dat=dat.s_mass[j];
      double m_par=pars[dat.n_distp+j];
      double z=(m_par-mean)/width;
      double ratio=pdf::s_norm(skew*z)/pdf::c_norm(skew*z);
    }
  }

  return 0;
}