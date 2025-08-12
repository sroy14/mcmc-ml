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

# Compilers
CXX     = g++
MPI_CXX = mpicxx

# Language & warnings
STD   = -std=c++20
WARN  = -Wall -Wextra -Wpedantic -Wshadow -Wformat=2 -Woverloaded-virtual -Wnon-virtual-dtor

# Python (for o2scl's <Python.h>)
PY_INCLUDES = -I/usr/include/python3.13 \
	-I/usr/lib/python3.13/site-packages/numpy/_core/include
PY_LDFLAGS  = -L/usr/lib  -ldl  -lm

# Build flags
CPPFLAGS = -Isrc $(PY_INCLUDES)
CXXFLAGS = $(STD) $(WARN) -O3 -DNDEBUG -pipe -MMD -MP
LDFLAGS  = $(PY_LDFLAGS)

# Layout
SRC_DIR = src
OBJDIR  = build

# Sources
SRCS    = $(wildcard $(SRC_DIR)/*.cpp)

# Objects (single build dir; distinct suffixes to avoid collisions)
OBJS_MPI   = $(patsubst $(SRC_DIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))
OBJS_NOMPI = $(patsubst $(SRC_DIR)/%.cpp,$(OBJDIR)/%.nompi.o,$(SRCS))

DEPS_MPI   = $(OBJS_MPI:.o=.d)
DEPS_NOMPI = $(OBJS_NOMPI:.o=.d)

# Executables
TARGET_MPI   = mc2ml
TARGET_NOMPI = mc2ml_nompi

# Default
all: help

# ===== MPI build =================================================
$(TARGET_MPI): $(OBJS_MPI)
	$(MPI_CXX) $(LDFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(MPI_CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# ===== Non-MPI build =============================================
$(TARGET_NOMPI): $(OBJS_NOMPI)
	$(CXX) $(LDFLAGS) -o $@ $^

$(OBJDIR)/%.nompi.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# ===== Convenience ===============================================
NP = 4

help:
	@echo "Targets:"
	@echo "  mc2ml         (MPI build, uses $(MPI_CXX))"
	@echo "  mc2ml_nompi   (non-MPI build, uses $(CXX))"
	@echo "  run           (mpirun -np \$$NP ./mc2ml)"
	@echo "  run_nompi     (./mc2ml_nompi)"
	@echo "  clean         (remove build/ and executables)"

run: $(TARGET_MPI)
	mpirun -np $(NP) ./$(TARGET_MPI)

run_nompi: $(TARGET_NOMPI)
	./$(TARGET_NOMPI)

clean:
	rm -rf $(OBJDIR) $(TARGET_MPI) $(TARGET_NOMPI)

# Auto deps
-include $(DEPS_MPI) $(DEPS_NOMPI)
