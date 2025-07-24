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
  along with Mcmc-ml. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/

#include "main.h"
#include <o2scl/hdf_io.h>
#include <o2scl/constants.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;
using namespace mc2ml;


int main::fill(const ubvector &pars, double wgt, 
               vector<double> &line, data &dat) {
  return 0;
}


int main::point(const ubvector &pars, std::ofstream &sout, 
                double &log_wgt, data &dat) {
  
  double mu=pars[0];
  double sigma=pars[1];
  double alpha=pars[2];
  
  double cf=1.0/sqrt(2.0*o2scl_const::pi)/sigma;

  log_wgt=0.0;

  for (double x=0.0; x<3.0; x+=0.1) {
    double pdf=exp(-0.5*pow((x-mu)/sigma, 2.0));
    double cdf=1.0+erf((x-mu)*alpha/sigma/sqrt(2.0));
    log_wgt+=log(cf*pdf*cdf);
  }

  return 0;

}

int main::deriv(const ubvector &pars, point_funct &pf, 
                ubvector &grad, data &dat, bool &success) {
  return 0;
}