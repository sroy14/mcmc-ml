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

#include "mcmc.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace mc2ml;


mcmc::mcmc() {

  mc_type="";
  ml_type="";
  
  set=make_shared<settings>();
  dat=make_shared<data>();
  
  m_ptr.resize(1);
  m_ptr[0]=new base;
  m_ptr[0]->set=set;
  m_ptr[0]->dat=dat;

} // mcmc()


int mcmc::set_threads(vector<string> &sv, bool itive_com) {
  
  if (sv.size()==1) {
    cerr << "Number of threads not specified in 'threads'." << endl;
    return 1;                                                          
  }
  if (mc_type.length()>0) {
    cerr << "Threads must be set before MC method." << endl;
    return 2;
  }
  if (ml_type.length()>0) {
    cerr << "Threads must be set before ML method." << endl;
    return 2;
  }

  size_t n_threads_old=n_threads;
  for(size_t i=0;i<n_threads_old;i++) {
    delete m_ptr[i];
  }

  n_threads=o2scl::stoszt(sv[1]);
  m_ptr.resize(n_threads);
  
  for (size_t i=0; i<n_threads; i++) {
    m_ptr[i]=new base;
    m_ptr[i]->set=set;
    m_ptr[i]->n_threads=n_threads;
  }
  
  return 0;

} // set_threads()


void mcmc::file_header(o2scl_hdf::hdf_file &hf) {

  mcmc_para_cli::file_header(hf);

  hf.set_szt("grid_size", set->grid_size);
  hf.sets("mc_method", mc_type);
  hf.sets("ml_method", ml_type);
  hf.seti("debug", set->debug);
  hf.setd("m_low", set->m_low);
  hf.setd("m_high", set->m_high);

  hdf_output(hf, dat->m_grid, "m_grid");

  return;

} // file_header()


int mcmc::mcmc_init() {

  /*cout << "(rank " << this->mpi_rank << ") "
       << "Start mcmc::mcmc_init()." << endl;*/
  
  if (m_ptr.size()<1) {
    O2SCL_ERR("mcmc::mcmc_init(): Object m_ptr invalid.", exc_esanity);
  }

  this->ret_value_counts.resize(this->n_threads);
  
  for(size_t it=0;it<this->n_threads;it++) {
    // The size must be at least (n_err_codes + 1)
    this->ret_value_counts[it].resize(4);
  }

  mcmc_para_cli::mcmc_init();

  // Add columns to table here
  
  if (set->inc_lmxb) {
    dat->load_data();
    for (size_t i=0; i<dat->n_stars; i++) {
      this->table->new_column("wgt_"+dat->s_names[i]);
    }
  } else {
    for (size_t i=0; i<set->grid_size; i++) {
      this->table->new_column("wgt_"+o2scl::itos(i));
    }
  }

  for (size_t i=0; i<n_threads; i++) {
    m_ptr[i]->dat->m_grid=o2scl::uniform_grid_end<double>
      (set->m_low, set->m_high, set->grid_size-1);
  }

  /*cout << "(rank " << this->mpi_rank << ") "
       << "End mcmc::mcmc_init()." << endl;*/

  return 0;

} // mcmc_init()


int mcmc::set_method_mc(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cerr << "MCMC method not specified." << endl;
    return o2scl::exc_efailed;
  }

  if (mc_type==sv[1]) {
    cerr << "MCMC method already set to " << sv[1] << "." << endl;
    return 0;
  }

  mc_type=sv[1];
  m_ptr[0]->mc_type=sv[1];

  return 0;

} // set_method_mc()


int mcmc::set_method_ml(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cerr << "ML method not specified." << endl;
    return o2scl::exc_efailed;
  }

  if (ml_type==sv[1]) {
    cerr << "ML method already set to " << sv[1] << "." << endl;
    return 0;
  }

  ml_type=sv[1];
  m_ptr[0]->ml_type=sv[1];

  return 0;

} // set_method_ml()


int mcmc::initial_point_best(vector<string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cerr << "Initial point file not specified." << endl;
    return 1;
  }

  size_t n_pars;
  string fname=sv[1];
  size_t pos=fname.find("<rank>");

  if (pos!=string::npos) {
    fname.replace(pos, 6, o2scl::itos(mpi_rank));
  }

  this->initial_points_file_best(fname, n_pars);

  return 0;

} // initial_points_file_best()


int mcmc::initial_point_last(vector<string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cerr << "Initial point file not specified." << endl;
    return 1;
  }

  size_t n_pars;
  string fname=sv[1];
  size_t pos=fname.find("<rank>");

  if (pos!=string::npos) {
    fname.replace(pos, 6, o2scl::itos(this->mpi_rank));
  }

  this->initial_points_file_last(fname, n_pars);

  return 0;

} // initial_points_file_last()


int mcmc::mcmc_func(vector<string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cerr << "MCMC method not specified." << endl;
    return 1;
  }

  vector<string> p_names, p_units;
  vector<string> d_names, d_units;

  vector<double> low, high, init;

  dat->get_param_info(p_names, p_units, low, high, set);
  set_names_units(p_names, p_units, d_names, d_units);

  if (this->initial_points.size()==0) {
    dat->set_init_point(init, set);
    this->initial_points.clear();
    ubvector ip(init.size());
    vector_copy(init, ip);
    this->initial_points.push_back(ip);
  }

  for(size_t i=0; i<p_names.size(); i++) {
    pvi.append(p_names[i]);
  }

  for(size_t i=0;i<m_ptr.size();i++) {
    m_ptr[i]->pvi=pvi;
  }

  size_t np=p_names.size();
  
  ubvector lo(low.size()), hi(high.size());
  vector_copy(low, lo);
  vector_copy(high, hi);

  vector<mc2ml::point_funct> pf(n_threads);
  vector<mc2ml::fill_funct>  ff(n_threads);

  using namespace std::placeholders;
  for (size_t i=0; i<n_threads; i++) {
    pf[i]=bind(mem_fn<int(const ubvector &, ofstream &, double &, data &)> 
              (&base::point), m_ptr[i], _2, ref(scr_out), _3, _4);
    ff[i]=bind(mem_fn<int(const ubvector &, double, vector<double> &, data &)>
              (&base::fill), m_ptr[i], _1, _2, _3, _4);
  }

  vector<data> dv(2*this->n_walk*this->n_threads);

  this->mcmc_fill(np, lo, hi, pf, ff, dv);

  return 0;

} // mcmc_func()


void mcmc::mcmc_setup_cli() {

  mcmc_para_cli::setup_cli(cl);
  set->setup_cli(cl);

  static const int n_opt=6;

  comm_option_s options[n_opt] = {
    /* Format: {short_opt, long_opt, help_desc, min_params, max_params,
     param_desc, help, func_ptr , com_type} */
    {'m', "mcmc", "Perform the MCMC simulation.", 0, 0, "",
      string("This is the main part of the code which performs ") +
       "the Markov Chain Monte Carlo simulation. Must set all " +
       "other options before this.",
      new comm_option_mfptr<mcmc>(this, &mcmc::mcmc_func),
      cli::comm_option_both},
    {'t', "threads", "Set number of OpenMP threads.", 1, 1, "<int>",
      string("The threads command must be set before ") +
       "any method selection and any changes to the settings.",
      new comm_option_mfptr<mcmc>(this, &mcmc::set_threads),
      cli::comm_option_both},
    {0, "method", "Set MCMC method.", 1, 1, "<string>",
      string("Choose the MCMC sampling method. Choices are: affine invariant ") +
       "('ai'), random walk ('rw'), Hamiltonian Monte Carlo ('hmc'), and " +
       "kernel density estimation ('kde').",
      new comm_option_mfptr<mcmc>(this, &mcmc::set_method_mc),
      cli::comm_option_both},
    {0, "model", "Set ML model.", 1, 1, "<string>",
      string("Choose the machine-learning model. ") +
       "Choices are: decision tree regressor ('dtr'), deep neural network " +
       "('dnn'), gaussian process ('gp'), and multi-layer perceptron " +
       "regressor ('mlpr').",
      new comm_option_mfptr<mcmc>(this, &mcmc::set_method_ml),
      cli::comm_option_both},
    {0, "initial-point-last", "Set initial point file name.", 1, 1, "<string>",
      string("Provide the name(s) of file(s) to read the initial points. Use ") +
       "this for continuing MCMC run using last points. If using MPI, format " +
       "the file name as, e.g., 'fname_<rank>'. ",
      new comm_option_mfptr<mcmc>(this, &mcmc::initial_point_last),
      cli::comm_option_both},
    {0, "initial-point-best", "Set initial point file name.", 1, 1, "<string>",
      string("Provide the name(s) of file(s) to read the initial points. Use ") +
       "this for restarting MCMC run using best points. If using MPI, format " +
       "the file name as, e.g., 'fname_<rank>'. ",
      new comm_option_mfptr<mcmc>(this,&mcmc::initial_point_best),
      cli::comm_option_both},
  };

  cl.set_comm_option_vec(n_opt, options);

  return;

} // setup_cli()