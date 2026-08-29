# SPDX-FileCopyrightText: 2026 Qiang Sun
# SPDX-License-Identifier: BSD-3-Clause

# Override FC/FFLAGS/LDLIBS on the command line when needed.
FC = gfortran
FFLAGS ?= -O2 -fopenmp -ffree-line-length-none
CHECKFLAGS ?= -O0 -fopenmp -ffree-line-length-none -fsyntax-only
BUILD_DIR ?= .build

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LDLIBS ?= -fopenmp -framework Accelerate
else
  LDLIBS ?= -fopenmp -llapack -lblas
endif

.PHONY: all syntax clean

all: au_ag_spectrum

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

syntax: | $(BUILD_DIR)
	$(FC) $(CHECKFLAGS) -J$(BUILD_DIR) Main_EM.f90

au_ag_spectrum: | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) Main_EM.f90 $(LDLIBS) -o $(BUILD_DIR)/$@

clean:
	rm -rf $(BUILD_DIR)
