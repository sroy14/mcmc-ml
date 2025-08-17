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

#include <o2scl/root_brent_gsl.h>

#include "data.h"

using namespace std;
using namespace o2scl;
using namespace mc2ml;


void data::get_param_info(vector<string> &names,
                          vector<string> &units,
                          vector<double> &low,
                          vector<double> &high,
                          shared_ptr<settings> set) {

  names.push_back("mean");
  units.push_back("Msun");
  low.push_back(0.5);
  high.push_back(3.0);

  names.push_back("log10_std");
  units.push_back("Msun");
  low.push_back(0.0);
  high.push_back(1.0);

  names.push_back("skewness");
  units.push_back("");
  low.push_back(-1.0);
  high.push_back(1.0);

  // LMXBs: n_stars = 5
  if (set->inc_lmxb) {
    for (size_t i=0; i<n_stars; i++) {
      names.push_back("M_"+s_names[i]);
      units.push_back(s_units[i]);
      low.push_back(0.5);
      high.push_back(2.5);
    }
  }

  return;

} // get_param_info()


void data::set_init_point(vector<double> &init,
                          shared_ptr<settings> set) {

  init.push_back(1.4);
  init.push_back(0.1);
  init.push_back(0.0);

  // LMXBs: n_stars = 5
  if (set->inc_lmxb) {
    for (size_t i=0; i<n_stars; i++) {
      init.push_back(s_mass[i]);
    }
  }

  return;

} // set_init_point()


void data::load_data() {

  // LMXBs: n_stars = 5
  vector<double> lo_68, hi_68;

  s_names.push_back("CygX2");
  s_mass.push_back(1.71);
  lo_68.push_back(0.21);

  s_names.push_back("XTEJ2123");
  s_mass.push_back(1.53);
  lo_68.push_back(0.42);

  s_names.push_back("4U1822");
  s_mass.push_back(1.96);
  lo_68.push_back(0.36);

  s_names.push_back("HerX1");
  s_mass.push_back(1.073);
  lo_68.push_back(0.36);

  s_names.push_back("2S0921");
  s_mass.push_back(1.44);
  lo_68.push_back(0.1);

  for (size_t i=0; i<n_stars; i++) {
    s_units.push_back("Msun");
  }

  solver s;
  for (size_t i=0; i<n_stars; i++) {
    hi_68.push_back(lo_68[i]); // Only for LMXBs
    double c=sqrt(hi_68[i]/lo_68[i]);
    double d=s.solve_1d(lo_68[i], hi_68[i]);
    c_68.push_back(c);
    d_68.push_back(d);
  }

}


double solver::funct_1d(double x, double &l, double &u) {
  double c=sqrt(u/l);
  return c*c*erf(u/(sqrt(2.0)*c*x)) - erf(-c*l/(sqrt(2.0)*x))
    - 0.68*(c*c+1.0);
}


double solver::solve_1d(double l, double u) {
  
  root_brent_gsl<> rbg;
  rbg.verbose=0;

  funct f=bind(mem_fn<double(double, double &, double &)>
		  (&solver::funct_1d), this, placeholders::_1, ref(l), ref(u));

  double x1=0.0, x2=1.0;

  rbg.solve_bkt(x1, x2, f);

  return x1;
}