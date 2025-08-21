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
  
  if (set->param_space=="S" || set->param_space=="M" ||
      set->param_space=="L" || set->param_space=="XL") { 
    // Include LMXBs

    names.push_back("mean_lmxb");
    units.push_back("Msun");
    low.push_back(0.5);
    high.push_back(2.5);

    names.push_back("log10std_lmxb");
    units.push_back("Msun");
    low.push_back(0.0);
    high.push_back(1.0);

    names.push_back("skewness_lmxb");
    units.push_back("");
    low.push_back(-1.0);
    high.push_back(1.0);

    for (int i=0; i<3; i++) {
      units.push_back("Msun");
      low.push_back(0.5);
      high.push_back(2.5);
    }

    // n_lmxb = 5
    for (size_t i=0; i<n_lmxb; i++) {
      names.push_back("M_"+s_names[i]);
      units.push_back(s_units[i]);
      low.push_back(0.5);
      high.push_back(2.5);
    }
  }

  if (set->param_space=="M" ||
      set->param_space=="L" ||
      set->param_space=="XL") { // Also include HMXBs

    names.push_back("mean_hmxb");
    units.push_back("Msun");
    low.push_back(0.5);
    high.push_back(2.5);

    names.push_back("log10std_hmxb");
    units.push_back("Msun");
    low.push_back(0.0);
    high.push_back(1.0);

    names.push_back("skewness_hmxb");
    units.push_back("");
    low.push_back(-1.0);
    high.push_back(1.0);

    // n_hmxb = 11
    for (size_t i=n_lmxb; i<n_lmxb+n_hmxb; i++) {
      names.push_back("M_"+s_names[i]);
      units.push_back(s_units[i]);
      low.push_back(0.5);
      high.push_back(2.5);
    }
  }

  if (set->param_space=="L" || set->param_space=="XL") {
    // Also include NS-NS
    names.push_back("mean_nsns");
    units.push_back("Msun");
    low.push_back(0.5);
    high.push_back(2.5);

    names.push_back("log10std_nsns");
    units.push_back("Msun");
    low.push_back(0.0);
    high.push_back(1.0);

    names.push_back("skewness_nsns");
    units.push_back("");
    low.push_back(-1.0);
    high.push_back(1.0);

    // n_nsns = 22
    for (size_t i=n_lmxb+n_hmxb; 
         i<n_lmxb+n_hmxb+n_nsns; i++) {
      names.push_back("M_"+s_names[i]);
      units.push_back(s_units[i]);
      low.push_back(0.5);
      high.push_back(2.5);
    }
  }

  if (set->param_space=="XL") { // Also include NS-WD

    names.push_back("mean_nswd");
    units.push_back("Msun");
    low.push_back(0.5);
    high.push_back(2.5);

    names.push_back("log10std_nswd");
    units.push_back("Msun");
    low.push_back(0.0);
    high.push_back(1.0);

    names.push_back("skewness_nswd");
    units.push_back("");
    low.push_back(-1.0);
    high.push_back(1.0);

    // n_nswd = 32
    for (size_t i=n_lmxb+n_hmxb+n_nsns; 
         i<n_lmxb+n_hmxb+n_nsns+n_nswd; i++) {
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
  
  if (set->param_space=="S" ||
      set->param_space=="M" ||
      set->param_space=="L" ||
      set->param_space=="XL") {
    init.push_back(1.4);
    init.push_back(-0.1);
    init.push_back(0.0);
    for (size_t i=0; i<n_lmxb; i++) {
      init.push_back(s_mass[i]);
    }
  }

  if (set->param_space=="M" ||
      set->param_space=="L" ||
      set->param_space=="XL") {
    init.push_back(1.4);
    init.push_back(-0.1);
    init.push_back(0.0);
    for (size_t i=n_lmxb; i<n_lmxb+n_hmxb; i++) {
      init.push_back(s_mass[i]);
    }
  }

  if (set->param_space=="L" ||
      set->param_space=="XL") {
    init.push_back(1.4);
    init.push_back(-0.1);
    init.push_back(0.0);
    for (size_t i=n_lmxb+n_hmxb; 
         i<n_lmxb+n_hmxb+n_nsns; i++) {
      init.push_back(s_mass[i]);
    }
  }

  if (set->param_space=="XL") {
    init.push_back(1.4);
    init.push_back(-0.1);
    init.push_back(0.0);
    for (size_t i=n_lmxb+n_hmxb+n_nsns; 
         i<n_lmxb+n_hmxb+n_nsns+n_nswd; i++) {
      init.push_back(s_mass[i]);
    }
  }

  return;

} // set_init_point()


void data::load_data(std::shared_ptr<settings> set) {

  vector<double> lo_68, hi_68;
  solver s;

  if (set->param_space=="S" || set->param_space=="M" ||
      set->param_space=="L" || set->param_space=="XL") {

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

    for (size_t i=0; i<n_lmxb; i++) {
      s_units.push_back("Msun");
      hi_68.push_back(lo_68[i]);
      double c=sqrt(hi_68[i]/lo_68[i]);
      double d=s.solve_1d(lo_68[i], hi_68[i]);
      c_68.push_back(c);
      d_68.push_back(d);
    }
  }

  if (set->param_space=="M" ||
      set->param_space=="L" ||
      set->param_space=="XL") {

    s_names.push_back("4U1700");
    s_mass.push_back(1.96); 
    lo_68.push_back(0.19); 

    s_names.push_back("SMCX1"); 
    s_mass.push_back(1.21); 
    lo_68.push_back(0.12); 

    s_names.push_back("CenX3"); 
    s_mass.push_back(1.57); 
    lo_68.push_back(0.16); 

    s_names.push_back("OAO1657"); 
    s_mass.push_back(1.74); 
    lo_68.push_back(0.3); 

    s_names.push_back("J013236");
    s_mass.push_back(2.0);
    lo_68.push_back(0.4);

    s_names.push_back("VelaX1");
    s_mass.push_back(2.12);
    lo_68.push_back(0.16);

    s_names.push_back("4U1538");
    s_mass.push_back(1.02);
    lo_68.push_back(0.17);

    s_names.push_back("LMCX4");
    s_mass.push_back(1.57);
    lo_68.push_back(0.11);

    s_names.push_back("EXO1722");
    s_mass.push_back(1.91);
    lo_68.push_back(0.45);

    s_names.push_back("SAXJ1802");
    s_mass.push_back(1.57);
    lo_68.push_back(0.25);

    s_names.push_back("XTEJ1855");
    s_mass.push_back(1.41);
    lo_68.push_back(0.24);

    for (size_t i=n_lmxb; i<n_lmxb+n_hmxb; i++) {
      s_units.push_back("Msun");
      hi_68.push_back(lo_68[i]);
      double c=sqrt(hi_68[i]/lo_68[i]);
      double d=s.solve_1d(lo_68[i], hi_68[i]);
      c_68.push_back(c);
      d_68.push_back(d);
    }
  }

  if (set->param_space=="L" ||
      set->param_space=="XL") {

    s_names.push_back("J0453p");
    s_mass.push_back(1.559);
    hi_68.push_back(0.004);
    lo_68.push_back(0.004);

    s_names.push_back("J0453c");
    s_mass.push_back(1.174);
    hi_68.push_back(0.004);
    lo_68.push_back(0.004);

    s_names.push_back("J1906p");
    s_mass.push_back(1.291);
    hi_68.push_back(0.011);
    lo_68.push_back(0.011);

    s_names.push_back("J1906c");
    s_mass.push_back(1.322);
    hi_68.push_back(0.011);
    lo_68.push_back(0.011);

    s_names.push_back("B1534p");
    s_mass.push_back(1.3332);
    hi_68.push_back(0.001);
    lo_68.push_back(0.001);

    s_names.push_back("B1534c");
    s_mass.push_back(1.3452);
    hi_68.push_back(0.001);
    lo_68.push_back(0.001);

    s_names.push_back("B1913p");
    s_mass.push_back(1.4398);
    hi_68.push_back(0.0002);
    lo_68.push_back(0.0002);

    s_names.push_back("B1913c");
    s_mass.push_back(1.3886);
    hi_68.push_back(0.0002);
    lo_68.push_back(0.0002);

    s_names.push_back("B2127p");
    s_mass.push_back(1.358);
    hi_68.push_back(0.01);
    lo_68.push_back(0.01);

    s_names.push_back("B2127c");
    s_mass.push_back(1.354);
    hi_68.push_back(0.01);
    lo_68.push_back(0.01);

    s_names.push_back("J0737A");
    s_mass.push_back(1.3381);
    hi_68.push_back(0.0007);
    lo_68.push_back(0.0007);

    s_names.push_back("J0737B");
    s_mass.push_back(1.2489);
    hi_68.push_back(0.0007);
    lo_68.push_back(0.0007);

    s_names.push_back("J1756p");
    s_mass.push_back(1.312);
    hi_68.push_back(0.017);
    lo_68.push_back(0.017);

    s_names.push_back("J1756c");
    s_mass.push_back(1.258);
    hi_68.push_back(0.017);
    lo_68.push_back(0.017);

    s_names.push_back("J1807p");
    s_mass.push_back(1.3655);
    hi_68.push_back(0.0021);
    lo_68.push_back(0.0021);

    s_names.push_back("J1807c");
    s_mass.push_back(1.2064);
    hi_68.push_back(0.002);
    lo_68.push_back(0.002);

    s_names.push_back("J1518p");
    s_mass.push_back(1.56);
    hi_68.push_back(0.13);
    lo_68.push_back(0.44);

    s_names.push_back("J1518c");
    s_mass.push_back(1.05);
    hi_68.push_back(0.45);
    lo_68.push_back(0.11);

    s_names.push_back("J1811p");
    s_mass.push_back(1.56);
    hi_68.push_back(0.24);
    lo_68.push_back(0.45);

    s_names.push_back("J1811c");
    s_mass.push_back(1.12);
    hi_68.push_back(0.47);
    lo_68.push_back(0.13);

    s_names.push_back("J1829p");
    s_mass.push_back(1.20);
    hi_68.push_back(0.12);
    lo_68.push_back(0.46);

    s_names.push_back("J1829c");
    s_mass.push_back(1.40);
    hi_68.push_back(0.46);
    lo_68.push_back(0.12);

    for (size_t i=n_lmxb+n_hmxb; 
        i<n_lmxb+n_hmxb+n_nsns; i++) {
      s_units.push_back("Msun");
      double c=sqrt(hi_68[i]/lo_68[i]);
      double d=s.solve_1d(lo_68[i], hi_68[i]);
      c_68.push_back(c);
      d_68.push_back(d);
    }
  }

  if (set->param_space=="XL") {

    s_names.push_back("J2045");
    s_mass.push_back(1.33);
    hi_68.push_back(0.3);
    lo_68.push_back(0.3);

    s_names.push_back("J2053");
    s_mass.push_back(1.40);
    hi_68.push_back(0.21);
    lo_68.push_back(0.21);

    s_names.push_back("J1713");
    s_mass.push_back(1.35);
    hi_68.push_back(0.07);
    lo_68.push_back(0.07);

    s_names.push_back("B1855");
    s_mass.push_back(1.37);
    hi_68.push_back(0.13); 
    lo_68.push_back(0.13);

    s_names.push_back("J0751");
    s_mass.push_back(1.72);
    hi_68.push_back(0.07);
    lo_68.push_back(0.07);

    s_names.push_back("J1141");
    s_mass.push_back(1.27);
    hi_68.push_back(0.01);
    lo_68.push_back(0.01);

    s_names.push_back("J1738");
    s_mass.push_back(1.47);
    hi_68.push_back(0.07);
    lo_68.push_back(0.07);

    s_names.push_back("J1614");
    s_mass.push_back(1.908);
    hi_68.push_back(0.016);
    lo_68.push_back(0.016);

    s_names.push_back("J0348");
    s_mass.push_back(2.01);
    hi_68.push_back(0.04);
    lo_68.push_back(0.04);

    s_names.push_back("J2222");
    s_mass.push_back(1.76);
    hi_68.push_back(0.06);
    lo_68.push_back(0.06);

    s_names.push_back("J2234");
    s_mass.push_back(1.393);
    hi_68.push_back(0.013);
    lo_68.push_back(0.013);

    s_names.push_back("J1949");
    s_mass.push_back(1.47);
    hi_68.push_back(0.43);
    lo_68.push_back(0.43);

    s_names.push_back("J1012");
    s_mass.push_back(1.83);
    hi_68.push_back(0.11);
    lo_68.push_back(0.11);

    s_names.push_back("J0437");
    s_mass.push_back(1.44);
    hi_68.push_back(0.07);
    lo_68.push_back(0.07);

    s_names.push_back("J1909");
    s_mass.push_back(1.48);
    hi_68.push_back(0.03);
    lo_68.push_back(0.03);

    s_names.push_back("J1802");
    s_mass.push_back(1.24);
    hi_68.push_back(0.11);
    lo_68.push_back(0.11);

    s_names.push_back("J1911");
    s_mass.push_back(1.34);
    hi_68.push_back(0.08);
    lo_68.push_back(0.08);

    s_names.push_back("J2043");
    s_mass.push_back(1.38);
    hi_68.push_back(0.13);
    lo_68.push_back(0.13);

    s_names.push_back("J0337");
    s_mass.push_back(1.4378);
    hi_68.push_back(0.0013);
    lo_68.push_back(0.0013);

    s_names.push_back("J1946");
    s_mass.push_back(1.828);
    hi_68.push_back(0.022);
    lo_68.push_back(0.022);

    s_names.push_back("J1918");
    s_mass.push_back(1.29);
    hi_68.push_back(0.1);
    lo_68.push_back(0.1);

    s_names.push_back("J1600");
    s_mass.push_back(2.3);
    hi_68.push_back(0.7);
    lo_68.push_back(0.7);

    s_names.push_back("J0621");
    s_mass.push_back(1.70);
    hi_68.push_back(0.10);
    lo_68.push_back(0.17);

    s_names.push_back("B2303");
    s_mass.push_back(1.38);
    hi_68.push_back(0.06);
    lo_68.push_back(0.10);

    s_names.push_back("J0024");
    s_mass.push_back(1.48);
    hi_68.push_back(0.03);
    lo_68.push_back(0.06);

    s_names.push_back("J0514");
    s_mass.push_back(1.49);
    hi_68.push_back(0.04);
    lo_68.push_back(0.27);

    s_names.push_back("B1516");
    s_mass.push_back(2.10);
    hi_68.push_back(0.19);
    lo_68.push_back(0.19);

    s_names.push_back("J1748I");
    s_mass.push_back(1.91);
    hi_68.push_back(0.02);
    lo_68.push_back(0.10);

    s_names.push_back("J1748J");
    s_mass.push_back(1.79);
    hi_68.push_back(0.02);
    lo_68.push_back(0.10);

    s_names.push_back("B1802");
    s_mass.push_back(1.26);
    hi_68.push_back(0.08);
    lo_68.push_back(0.17);

    s_names.push_back("B1911");
    s_mass.push_back(1.40);
    hi_68.push_back(0.16);
    lo_68.push_back(0.10);

    s_names.push_back("J0740");
    s_mass.push_back(2.14);
    hi_68.push_back(0.20);
    lo_68.push_back(0.18);

    for (size_t i=n_lmxb+n_hmxb+n_nsns; 
         i<n_lmxb+n_hmxb+n_nsns+n_nswd; i++) {
      s_units.push_back("Msun");
      double c=sqrt(hi_68[i]/lo_68[i]);
      double d=s.solve_1d(lo_68[i], hi_68[i]);
      c_68.push_back(c);
      d_68.push_back(d);
    }
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