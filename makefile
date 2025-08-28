#  -------------------------------------------------------------------
#  
#  Copyright (C) 2012-2025, Mahmudul Hasan Anik, Satyajit Roy, and 
#  Andrew W. Steiner
#  
#  This file is part of mcmc-ml.
#
#  Mcmc-ml is free software; you can redistribute it and/or modify it 
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Mcmc-ml is distributed in the hope that it will be useful, but 
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with mcmc-ml. If not, see <http://www.gnu.org/licenses/>.
#
#  -------------------------------------------------------------------

# ---------------------------------------------------------------
# mc2ml Makefile (src -> build)
# Build one at a time: mc2ml (non-MPI) | mc2ml_mpi (MPI)
# ---------------------------------------------------------------

# Default C++ compilers
CXX = g++
CXX_MPI = mpic++

# C++ Standard
STD = -std=c++20

# Warning flags
WARN = -Wall -Wno-unused -Wshadow -Wformat=2 -Woverloaded-virtual

# Optimization flags
OPT_FLAGS = -O3 -DNDEBUG -pipe -MMD -MP

# OpenMP flags
OMP_FLAGS = -fopenmp -DO2SCL_OPENMP

# MPI flags
MPI_FLAGS = -DO2SCL_MPI

# HDF5 compression & header
HDF_FLAGS = -DO2SCL_HDF5_COMP -DO2SCL_PLAIN_HDF5_HEADER

# Python includes and libraries
PY_INC = -I/usr/include/python3.13 \
	-I/usr/lib/python3.13/site-packages/numpy/_core/include
PY_LIB = -L/usr/lib -ldl -lm

# Compiler flags
CFLAGS = $(OMP_FLAGS) $(HDF_FLAGS) $(PY_INC)
CFLAGS_MPI = $(CFLAGS) $(MPI_FLAGS)
CXXFLAGS = $(STD) $(OPT_FLAGS) $(WARN)

# Library flags
LIB_FLAGS = -lo2scl -lhdf5 -lgsl -lreadline -lpython3.13
LDFLAGS = -fopenmp $(PY_LIB) $(LIB_FLAGS)

# Layout
SRC_DIR = src
OBJ_DIR = build
OUT_DIR = out

# Sources
SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# Objects & deps (single build dir; MPI uses .mpi.o/.mpi.d)
OBJS      = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
OBJS_MPI  = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.mpi.o,$(SRCS))
DEPS      = $(OBJS:.o=.d)
DEPS_MPI  = $(OBJS_MPI:.o=.d)   # -> .mpi.d

# Executables
TARGET     = mc2ml
TARGET_MPI = mc2ml_mpi

# Default
all: help

$(OUT_DIR):
	@mkdir -p $(OUT_DIR)

# ----- Non-MPI build ---------------------------------------------
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -c $< -o $@

# ----- MPI build -------------------------------------------------
$(TARGET_MPI): $(OBJS_MPI)
	$(CXX_MPI) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.mpi.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX_MPI) $(CXXFLAGS) $(CFLAGS_MPI) -c $< -o $@

# ----- Convenience ----------------------------------------------
NP = 8

help:
	@echo "Targets:"
	@echo "  $(TARGET)       (non-MPI)"
	@echo "  $(TARGET_MPI)   (MPI)"
	@echo "  run             (./$(TARGET))"
	@echo "  run_MPI         (mpirun -np \$$NP ./$(TARGET_MPI))"
	@echo "  info            (print -L/-l entries in flags)"
	@echo "  clean           (remove build/ and executables)"

run: $(TARGET) $(OUT_DIR)
	./$(TARGET) -threads 1 \
	-method hmc -param-space S -set prefix $(OUT_DIR)/hmc -set n_walk 24 \
	-set max_iters 10000 -set file_update_time 1 \
	-set verbose 1 -set mcmc_verbose 2  \
	-mcmc

run_mpi: $(TARGET_MPI) $(OUT_DIR)
	mpirun -np $(NP) ./$(TARGET_MPI) -threads 1 \
	-method hmc -param-space S -set prefix $(OUT_DIR)/hmc -set n_walk 24 \
	-set max_iters 10000 -set file_update_time 1 \
	-set verbose 1 -set mcmc_verbose 2  \
	-mcmc

# Print library-related tokens from flags
info:
	@printf "`tput bold`# Default C++ compilers\n`tput sgr0`"
	@echo "CXX = $(CXX)"
	@echo "CXX_MPI = $(CXX_MPI)"
	@echo

	@printf "`tput bold`# C++ Standard\n`tput sgr0`"
	@echo "STD = $(STD)"
	@echo

	@printf "`tput bold`# Warning flags\n`tput sgr0`"
	@echo "WARN = $(WARN)"
	@echo

	@printf "`tput bold`# Optimization flags\n`tput sgr0`"
	@echo "OPT_FLAGS = $(OPT_FLAGS)"
	@echo

	@printf "`tput bold`# OpenMP/MPI flags\n`tput sgr0`"
	@echo "OMP_FLAGS = $(OMP_FLAGS)"
	@echo "MPI_FLAGS = $(MPI_FLAGS)"
	@echo

	@printf "`tput bold`# HDF5 compression & header flags\n`tput sgr0`"
	@echo "HDF_FLAGS = $(HDF_FLAGS)"
	@echo

	@printf "`tput bold`# Python includes and libraries\n`tput sgr0`"
	@echo "PY_INC = $(PY_INC)"
	@echo "PY_LIB = $(PY_LIB)"
	@echo

	@printf "`tput bold`# Compiler flags\n`tput sgr0`"
	@echo "CXXFLAGS = $(STD) $(OPT_FLAGS) $(WARN)"
	@echo "	   $(OMP_FLAGS) $(HDF_FLAGS) $(MPI_FLAGS)"
	@echo "           $(PY_INC)"
	@echo

	@printf "`tput bold`# Library flags\n`tput sgr0`"
	@echo "LDFLAGS = -fopenmp $(PY_LIB) $(LIB_FLAGS)"
	@echo

	@printf "`tput bold`# Directory layout\n`tput sgr0`"
	@echo "SRC_DIR = $(SRC_DIR)/"
	@echo "OBJ_DIR = $(OBJ_DIR)/"
	@echo "OUT_DIR = $(OUT_DIR)/"
	@echo

	@printf "`tput bold`# Executables (MPI/non-MPI)\n`tput sgr0`"
	@echo "TARGETS = $(TARGET_MPI) $(TARGET)"

clean:
	rm -rf $(OBJ_DIR) $(TARGET) $(TARGET_MPI)

# Auto deps
-include $(DEPS) $(DEPS_MPI)
