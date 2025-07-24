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

#ifndef DATA_H
#define DATA_H

#include <o2scl/table.h>

namespace mc2ml {

  class data {

  public:

    o2scl::table<> grid;

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

  };

}

#endif // DATA_H