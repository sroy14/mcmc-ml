help:
	@echo "Makefile targets:"
	@echo "─────────────────────────────────────────────────────────"
	@echo "ex_mcmc:          later."
	@echo "ex_mcmc_kde:      later."
	@echo "ex_mcmc_new:      later."
	@echo "ex_mcmc_nn:       later."

# ----------------------------------------------------------------
# Various user-specific settings
# ----------------------------------------------------------------
# LIBS is the list of libraries
# LCXX is the local C++ compiler
# LCFLAGS are the local C++ compiler flags

# Default settings
# LCXX = $(CXX)
# LMPI_CXX = $(MPI_CXX)
# LIBS = -lo2scl -lhdf5 -lgsl -lreadline $(LDFLAGS) 
# LMPI_CFLAGS = -O3 -DO2SCL_MPI -fopenmp $(CFLAGS) $(MPI_CFLAGS)
# LCFLAGS = -O3 -std=c++11 -DNO_MPI -fopenmp $(CFLAGS)

# ----------------------------------------------------------------
# UTK-specific settings
# ----------------------------------------------------------------

ifdef UTKNA_MAKEFILE

include $(UTKNA_MAKEFILE)

# UTK configuration

LIBS = $(UTKNA_O2SCL_LIBS)
INCS = $(UTKNA_O2SCL_INCS)
LCXX = $(UTKNA_CXX)
LMPI_CXX = $(UTKNA_MPI_CXX)
LCFLAGS = $(UTKNA_CFLAGS) $(INCS)
LMPI_CFLAGS = $(UTKNA_MPI_CFLAGS) $(UTKNA_OPENMP_FLAGS) $(INCS)

endif

# ----------------------------------------------------------------
# Main
# ----------------------------------------------------------------
# Object files - openmp no mpi
# ----------------------------------------------------------------
ex_mcmc.o: ex_mcmc.cpp  
	$(LCXX) $(LCFLAGS) -o ex_mcmc.o -c ex_mcmc.cpp 

ex_mcmc_kde.o: ex_mcmc_kde.cpp  
	$(LCXX) $(LCFLAGS) -o ex_mcmc_kde.o -c ex_mcmc_kde.cpp 

ex_mcmc_new.o: ex_mcmc_new.cpp  
	$(LCXX) $(LCFLAGS) -o ex_mcmc_new.o -c ex_mcmc_new.cpp 

ex_mcmc_nn.o: ex_mcmc_nn.cpp  
	$(LCXX) $(LCFLAGS) -o ex_mcmc_nn.o -c ex_mcmc_nn.cpp 

# ----------------------------------------------------------------
# Executables - openmp no mpi
# ----------------------------------------------------------------

ex_mcmc: ex_mcmc.o 
	$(LCXX) $(LCFLAGS) -o ex_mcmc ex_mcmc.o $(LIBS)

ex_mcmc_kde: ex_mcmc_kde.o 
	$(LCXX) $(LCFLAGS) -o ex_mcmc_kde ex_mcmc_kde.o $(LIBS)

ex_mcmc_new: ex_mcmc_new.o 
	$(LCXX) $(LCFLAGS) -o ex_mcmc_new ex_mcmc_new.o $(LIBS)

ex_mcmc_nn: ex_mcmc_nn.o 
	$(LCXX) $(LCFLAGS) -o ex_mcmc_nn ex_mcmc_nn.o $(LIBS)

clean:
	rm -f *.o *_out *_scr ex_mcmc ex_mcmc_kde ex_mcmc_new ex_mcmc_nn 

# This optional file, makefile.user, is an alternate place to store
# the user's makefile targets.
#-include makefile.sroy