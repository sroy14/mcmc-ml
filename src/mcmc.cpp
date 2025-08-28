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

#include <iostream>
#include <typeinfo>

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
  
  for (size_t i=0; i<dat->n_stars; i++) {
    this->table->new_column("wgt_"+dat->s_names[i]);
  }

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

  /*if (mc_type=="ai") this->aff_inv=1;
  else this->aff_inv=0;*/

  return 0;

} // set_method_mc()


int mcmc::set_method_ml(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cout << "ML method not specified: running standalone MCMC." 
         << endl;
    return 0;
  }

  if (ml_type==sv[1]) {
    cerr << "ML method already set to " << sv[1] << "." << endl;
    return 0;
  }

  ml_type=sv[1];
  m_ptr[0]->ml_type=sv[1];

  return 0;

} // set_method_ml()


int mcmc::set_param_space(vector<string> &sv, bool itive_com) {

  cout << "mcmc::set_param_space()" << endl;

  if (sv.size()<2) {
    cerr << "Parameter space not specified." << endl;
    return o2scl::exc_einval;
  }

  if (sv[1]=="S") {
    set->n_pops=1;
    n_params=3+dat->n_lmxb;
  }
  else if (sv[1]=="M") {
    set->n_pops=2;
    n_params=6+dat->n_lmxb+dat->n_hmxb;
  }
  else if (sv[1]=="L") {
    set->n_pops=3;
    n_params=9+dat->n_lmxb+dat->n_hmxb+dat->n_nsns;
  }
  else if (sv[1]=="XL") {
    set->n_pops=4;
    n_params=12+dat->n_lmxb+dat->n_hmxb+dat->n_nsns+dat->n_nswd;
  }
  else {
    cerr << "Invalid parameter space: " << sv[1] << endl;
    return o2scl::exc_einval;
  }

  return 0;

} // set_param_space()


int mcmc::initial_point_best(vector<string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cerr << "Initial point file not specified." << endl;
    return o2scl::exc_einval;
  }

  string fname=sv[1];
  size_t pos=fname.find("<rank>");

  if (pos!=string::npos) {
    fname.replace(pos, 6, o2scl::itos(mpi_rank));
  }

  this->initial_points_file_best(fname, n_params);

  return 0;

} // initial_points_file_best()


int mcmc::initial_point_last(vector<string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cerr << "Initial point file not specified." << endl;
    return o2scl::exc_einval;
  }

  string fname=sv[1];
  size_t pos=fname.find("<rank>");

  if (pos!=string::npos) {
    fname.replace(pos, 6, o2scl::itos(this->mpi_rank));
  }

  this->initial_points_file_last(fname, n_params);

  return 0;

} // initial_points_file_last()


int mcmc::dump_params(std::vector<std::string>&, bool) {
  std::cout << "Parameters registered (" << cl.par_list.size() << ")\n";
  for (const auto &kv : cl.par_list) {
    const char *ty = "unknown";
    if (dynamic_cast<o2scl::cli::parameter_size_t*>(kv.second))   ty = "size_t";
    else if (dynamic_cast<o2scl::cli::parameter_double*>(kv.second))  ty = "double";
    else if (dynamic_cast<o2scl::cli::parameter_int*>(kv.second))     ty = "int";
    else if (dynamic_cast<o2scl::cli::parameter_bool*>(kv.second))    ty = "bool";
    else if (dynamic_cast<o2scl::cli::parameter_string*>(kv.second))  ty = "string";
    std::cout << "  " << kv.first << " : " << ty << "\n";
  }
  return 0;
}


int mcmc::print_config(std::vector<std::string>&, bool) {
  std::cout
    << "CONFIG:\n"
    << "  prefix=" << prefix << "\n"
    << "  n_walk=" << n_walk << "\n"
    << "  max_iters=" << max_iters << "\n"
    << "  file_update_time=" << file_update_time << "\n"
    << "  verbose=" << set->verbose << "\n"
    << "  mcmc_verbose=" << this->verbose << "\n";
  return 0;
}


int mcmc::mcmc_func(vector<string> &sv, bool itive_com) {

  if (mc_type.empty()) {
    cerr << "MCMC method not specified." << endl;
    return 1;
  }
  
  // Initialize grid before passing 'dat' object to base::point()
  dat->m_grid=o2scl::uniform_grid_end<double>
      (set->m_low, set->m_high, set->grid_size-1);

  vector<string> p_names, p_units;
  vector<string> d_names, d_units;

  vector<double> low, high, init;

  dat->load_data(set);
  dat->get_param_info(p_names, p_units, low, high, set);
  set_names_units(p_names, p_units, d_names, d_units);

  if (this->initial_points.size()==0) {
    dat->set_init_point(init, set);
    this->initial_points.clear();
    ubvector ip(init.size());
    vector_copy(init, ip);
    this->initial_points.push_back(ip);
  }

  for(size_t i=0; i<n_params; i++) {
    pvi.append(p_names[i]);
  }

  for(size_t i=0;i<m_ptr.size();i++) {
    m_ptr[i]->pvi=pvi;
  }

  size_t np=n_params;
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

  size_t sz=2*this->n_walk*this->n_threads;
  vector<data> dv(sz, *dat);

  if (mc_type=="hmc") {

#ifdef O2SCL_MPI
    if (mpi_size>1 && mpi_rank>=1) {
      int tag=0, buffer=0;
      MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,
         tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
#endif

    shared_ptr<mcmc_stepper_hmc<point_funct,data,ubvector>> hmc_stepper
      (new mcmc_stepper_hmc<point_funct,data,ubvector>);
    stepper=hmc_stepper;

    hmc_stepper->allocate(n_threads, n_params);

    hmc_stepper->auto_grad.resize(np);
    for (size_t i=0; i<np; i++) {
      hmc_stepper->auto_grad[i]=false;
    }

    hmc_stepper->hmc_step.resize(np);
    for (size_t i=0; i<np; i++) {
      hmc_stepper->hmc_step[i]=1.0e-2*(high[i]-low[i]);
    }

    hmc_stepper->traj_length=1;

    vector<mc2ml::deriv_funct> gf(n_threads);
    for (size_t i=0; i<n_threads; i++) {
      gf[i]=bind(mem_fn<int(const ubvector &, point_funct &, ubvector &, 
                            data &)>(&base::deriv), m_ptr[i], _2, _3, _4, _5);
    }

    hmc_stepper->set_gradients(gf);

#ifdef O2SCL_MPI
    if (mpi_size>1 && mpi_rank>=1) {
      int tag=0, buffer=0;
      MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
         tag,MPI_COMM_WORLD);
    }
#endif

  }

  this->mcmc_fill(np, lo, hi, pf, ff, dv);

  return 0;

} // mcmc_func()


void mcmc::mcmc_setup_cli() {

  mcmc_para_cli::setup_cli(cl);
  set->setup_cli(cl);

  static const int n_opt=9;

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
    {'p', "param-space", "Set parameter space.", 1, 1, "<string>",
      string("Specify the size of the parameter space for MCMC. ") +
       "Choices are: 'S', 'M', 'L', 'XL'.",
      new comm_option_mfptr<mcmc>(this, &mcmc::set_param_space),
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
    {0, "dump-params", "List all -set parameter names and types.", 0, 0, "", "",
      new comm_option_mfptr<mcmc>(this, &mcmc::dump_params),
      cli::comm_option_both},
    {0, "print-config", "Show current run configuration.", 0, 0, "", "",
      new comm_option_mfptr<mcmc>(this, &mcmc::print_config),
      cli::comm_option_both},

  };

  cl.set_comm_option_vec(n_opt, options);

  return;

} // setup_cli()