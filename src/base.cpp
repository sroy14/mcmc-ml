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
  
  for (size_t i=0; i<dat.n_stars; i++) {
    line.push_back(dat.wgt_star.get("wgt", i));
  }

  return 0;

}


int base::point(const ubvector &pars, std::ofstream &sout, 
                double &log_wgt, data &dat) {
  
  log_wgt=0.0;

  if (dat.wgt_star.get_ncolumns()==0) {
    dat.wgt_star.new_column("wgt");
    dat.wgt_star.set_nlines(dat.n_stars);
  }

  vector<size_t> v={dat.n_lmxb, dat.n_hmxb, dat.n_nsns, dat.n_nswd};

  for (size_t i=0,j=0; i<set->n_pops; j+=v[i],i++) {
    double mean=pars[3*i+0];
    double width=pow(10.0, pars[3*i+1]);
    double skew=pars[3*i+2];

    for (size_t k=0; k<v[i]; k++) {
      double m_dat=dat.s_mass[j+k];
      double asym=dat.c_68[j+k];
      double scale=dat.d_68[j+k];
      double m_par=pars[3*set->n_pops+j+k];
      double wgt=pdf::asym_norm(m_dat-m_par, asym, scale) * 
                pdf::skewed_norm(m_par, mean, width, skew);

      if (wgt<=0.0) {
        sout << "base::point(): Star " << dat.s_names[j+k] 
              << " returned zero weight." << endl;
        // log_wgt=double(dat.ix_wgt_zero)-100.0;
        return dat.ix_wgt_zero;
      }

      dat.wgt_star.set("wgt", j+k, wgt);
      log_wgt+=log(wgt);
    }
  }

  return 0;

}


int base::deriv(const ubvector &pars, point_funct &pf, 
                ubvector &grad, data &dat) {

  double log_wgt=0.0;
  size_t np=pars.size();
  int f_ret=pf(np, pars, log_wgt, dat);
  if (f_ret!=0) return dat.ix_grad_failure;

  vector<size_t> v={dat.n_lmxb, dat.n_hmxb, dat.n_nsns, dat.n_nswd};

  if (grad.size()!=np) grad.resize(np);

  for (size_t i=0,j=0; i<set->n_pops; j+=v[i],i++) {
    double mean=pars[3*i+0];
    double width=pow(10.0, pars[3*i+1]);
    double skew=pars[3*i+2];

    double sum_mean=0.0;
    double sum_width=0.0;
    double sum_skew=0.0;

    for (size_t k=0; k<v[i]; k++) {
      double m_dat=dat.s_mass[j+k];
      double asym=dat.c_68[j+k];
      double scale=dat.d_68[j+k];
      double m_par=pars[3*set->n_pops+j+k];

      double expo=(m_par-mean)/width;
      double ncdf=pdf::c_norm(skew*expo);
      double npdf=pdf::s_norm(skew*expo);
      double ratio=npdf/ncdf;

      double d_mean=(1.0/width)*(expo-skew*ratio);
      double d_width=(1.0/width)*((expo*expo-1.0)-skew*expo*ratio);
      double d_skew=(expo*ratio);
      double dm_sn=(1.0/width)*(-expo+skew*ratio);

      sum_mean+=d_mean;
      sum_width+=d_width;
      sum_skew+=d_skew;

      double f_an=pdf::asym_norm(m_dat-m_par, asym, scale);
      double dm_an=(-pdf::dan_dx(m_dat-m_par, asym, scale))/f_an;

      grad[3*set->n_pops+j+k]=-(dm_sn+dm_an);
    }

    grad[3*i+0]=-sum_mean;
    grad[3*i+1]=-(log(10.0)*width)*sum_width;
    grad[3*i+2]=-sum_skew;
  }

  return 0;
}