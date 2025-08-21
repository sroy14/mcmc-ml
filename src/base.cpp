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

  if (set->param_space=="S" || set->param_space=="M" ||
      set->param_space=="L" || set->param_space=="XL") {
    for (size_t i=0; i<dat.n_lmxb; i++) {
      line.push_back(dat.wgt_star.get("wgt", i));
    }
  }

  if (set->param_space=="M" ||
      set->param_space=="L" ||
      set->param_space=="XL") {
    for (size_t i=dat.n_lmxb; i<dat.n_lmxb+dat.n_hmxb; i++) {
      line.push_back(dat.wgt_star.get("wgt", i));
    }
  }

  if (set->param_space=="L" || set->param_space=="XL") {
    for (size_t i=dat.n_lmxb+dat.n_hmxb; 
         i<dat.n_lmxb+dat.n_hmxb+dat.n_nsns; i++) {
      line.push_back(dat.wgt_star.get("wgt", i));
    }
  }

  if (set->param_space=="XL") {
    for (size_t i=dat.n_lmxb+dat.n_hmxb+dat.n_nsns; 
         i<dat.n_lmxb+dat.n_hmxb+dat.n_nsns+dat.n_nswd; i++) {
      line.push_back(dat.wgt_star.get("wgt", i));
    }
  }

  return 0;

}


int base::point(const ubvector &pars, std::ofstream &sout, 
                double &log_wgt, data &dat) {
  
  log_wgt=0.0;
  size_t n_stars;

  if (set->param_space=="S") {
    n_stars=dat.n_lmxb;
  } else if (set->param_space=="M") {
    n_stars=dat.n_lmxb+dat.n_hmxb;
  } else if (set->param_space=="L") {
    n_stars=dat.n_lmxb+dat.n_hmxb+dat.n_nsns;
  } else if (set->param_space=="XL") {
    n_stars=dat.n_lmxb+dat.n_hmxb+dat.n_nsns+dat.n_nswd;
  } else {
    sout << "base::point(): Unknown parameter space." << endl;
    return -1;
  }

  vector<double> mean(n_stars), width(n_stars);
  vector<double> skew(n_stars), mass(n_stars);

  if (dat.wgt_star.get_ncolumns()==0) {
    dat.wgt_star.new_column("wgt");
    dat.wgt_star.set_nlines(n_stars);
  }

  if (set->inc_lmxb) {

    double mean=pars[pvi["mean"]];
    double width=pow(10.0, pars[pvi["log10_std"]]);
    double skew=pars[pvi["skewness"]];

    for (size_t i=0; i<n_stars; i++) {

      double m_dat=dat.s_mass[i];
      double asym=dat.c_68[i];
      double scale=dat.d_68[i];
      double m_par=pars[3+i];
      double wgt=pdf::asym_norm(m_dat-m_par, asym, scale) * 
                 pdf::skewed_norm(m_par, mean, width, skew);

      if (wgt<=0.0) {
        sout << "base::point(): Star " << dat.s_names[i] 
             << " returned zero weight." << endl;
        // log_wgt=double(dat.ix_wgt_zero)-100.0;
        return dat.ix_wgt_zero;
      }

      dat.wgt_star.set("wgt", i, wgt);
      log_wgt+=log(wgt);

    }
  }

  return 0;

}

int base::deriv(const ubvector &pars, point_funct &pf, 
                ubvector &grad, data &dat) {
  
  /*size_t n_params=pars.size();
  if (grad.size()!=n_params) grad.resize(n_params);

  double weights=0.0;

  for (size_t i=0; i<dat.n_stars; i++) {
    weights+=dat.wgt_star.get("wgt", i);
  }

  vector<double> sum_m(dat.n_pops,0.0);
  vector<double> sum_s(dat.n_pops,0.0);
  vector<double> sum_a(dat.n_pops,0.0);

  // Main loop over all (i,j)
  for (size_t i=0; i<dat.n_pops; i++) {
    
    double mean=pars[pvi["mean"]];
    double width=pow(10.0, pars[pvi["log10_std"]]);
    double skew=pars[pvi["skewness"]];

    for (size_t j=0; j<dat.n_lmxb; j++) {
      double m_dat=dat.s_mass[j];
      double m_par=pars[dat.n_distp+j];
      double expo=(m_par-mean)/width;
      double ratio=pdf::s_norm(skew*expo)/pdf::c_norm(skew*expo);
      double ddm_sn=(1.0/width)*(expo-skew*ratio);
      double ddw_sn=(1.0/width)*((expo*expo-1.0)-skew*expo*ratio);
      double dds_sn=expo*ratio;
      double ddM_sn=(1.0/width)*(-expo+skew*ratio);
      sum_m[i]+=ddm_sn;
      sum_s[i]+=ddw_sn;
      sum_a[i]+=dds_sn;
      double w=dat.s_mass[i]-m_par;
      double c=dat.c_68[i];
      double d=dat.d_68[i];
      double f_an=pdf::asym_norm(w,c,d);
      double ddM_an=(-pdf::dan_dx(w,c,d))/f_an;
      grad[i]=weights*(ddM_sn+ddM_an);
    }
  }*/

  return 0;
}