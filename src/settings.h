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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <o2scl/cli.h>

namespace mc2ml {

  class settings {

  public:

    settings() {
      grid_size=100;
      verbose=0;
      debug=false;
      m_low=0.2;
      m_high=3.0;
      data_dir="data";
      n_pops=1;
    }

    o2scl::cli::parameter_size_t p_grid_size;
    o2scl::cli::parameter_int p_verbose;
    o2scl::cli::parameter_bool p_debug;
    o2scl::cli::parameter_double p_m_low;
    o2scl::cli::parameter_double p_m_high;
    o2scl::cli::parameter_string p_data_dir;

    size_t grid_size;
    size_t n_pops;
    int verbose;
    bool debug;
    double m_low;
    double m_high;
    std::string data_dir;

    void setup_cli(o2scl::cli &cl) {
      p_grid_size.s=&grid_size;
      p_grid_size.help="Grid size (default 100).";
      cl.par_list.insert(std::make_pair("grid_size",&p_grid_size));

      p_verbose.i=&verbose;
      p_verbose.help="Controls verbosity (default 0).";
      cl.par_list.insert(std::make_pair("verbose",&p_verbose));

      p_debug.b=&debug;
      p_debug.help="If true, output debugging info (default false).";
      cl.par_list.insert(std::make_pair("debug",&p_debug));

      p_m_low.d=&m_low;
      p_m_low.help="Smallest mass grid point in Msun (default 0.2).";
      cl.par_list.insert(std::make_pair("m_low",&p_m_low));
      
      p_m_high.d=&m_high;
      p_m_high.help="Largest mass grid point in Msun (default 3.0).";
      cl.par_list.insert(std::make_pair("m_high",&p_m_high));

      p_data_dir.str=&data_dir;
      p_data_dir.help="Directory for input data files (default 'data').";
      cl.par_list.insert(std::make_pair("data_dir",&p_data_dir));

      return;
    }

  }; // class settings
} // namespace mc2ml

#endif // SETTINGS_H