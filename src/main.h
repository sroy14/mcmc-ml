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

#ifndef MAIN_H
#define MAIN_H

#include "data.h"

namespace mc2ml {

  typedef boost::numeric::ublas::vector<double> ubvector;
  
  typedef std::function<int(size_t, const ubvector &, 
                            double &, data &)> point_funct;
  
  typedef std::function<int(const ubvector &, double, 
                            std::vector<double> &, data &)> fill_funct;

  typedef std::function<int(size_t, const ubvector &, point_funct &, 
                            ubvector &, data &, bool &)> deriv_funct;
  
  class main {
  
  public:
    
    int n_threads;
    
    std::string mcmc_method;

    main() {
      n_threads=1;
    }

    virtual ~main() {
    }

    virtual int point(const ubvector &, std::ofstream &, double &, data &);
    virtual int deriv(const ubvector &, point_funct &, ubvector &, data &, bool &);
    virtual int fill(const ubvector &, double, std::vector<double> &, data &);

  };


}

#endif // MAIN_H