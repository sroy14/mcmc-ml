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

#ifndef MCMC_H
#define MCMC_H

#ifdef O2SCL_MPI
#include <mpi.h>
#endif

#ifdef O2SCL_PYTHON
#include <Python.h>
#endif

#include <o2scl/mcmc_para.h>
#include <o2scl/hdf_file.h>
#include <o2scl/kde_python.h>
#include <o2scl/cli.h>

#include "base.h"

namespace mc2ml {
  
  class mcmc :
    public o2scl::mcmc_para_cli
      <point_funct, fill_funct, data, ubvector> {

  protected:
    
    std::string mc_type;
    std::string ml_type;
    std::vector<base*> m_ptr;
    std::shared_ptr<settings> set;
    std::shared_ptr<data> dat;
    o2scl::vec_index pvi;
    size_t n_params;

    virtual int set_method_mc(std::vector<std::string> &, bool);
    virtual int set_method_ml(std::vector<std::string> &, bool);
    virtual void file_header(o2scl_hdf::hdf_file &);
    virtual int mcmc_init();
    virtual int mcmc_func(std::vector<std::string> &, bool);
    virtual int set_threads(std::vector<std::string> &, bool);
    virtual int set_param_space(std::vector<std::string> &, bool);
    virtual int initial_point_last(std::vector<std::string> &, bool);
    virtual int initial_point_best(std::vector<std::string> &, bool);

  public:

    o2scl::cli cl;

    mcmc();

    virtual ~mcmc() {}

    virtual void mcmc_setup_cli();

  }; // class mcmc

} // namespace mc2ml


#endif // MCMC_H