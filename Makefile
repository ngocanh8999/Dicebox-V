# ==============================================================================
# Makefile for DICEBOX-V
# Computer Physics Communications (CPC) Program Library
#
# ==============================================================================
# USAGE INSTRUCTIONS:
# ------------------------------------------------------------------------------
# 1. make / make all   : Build ALL binaries (DICEBOX-V OpenMP & DICEBOX-1.0.1)
#                        -> Produces: bin/dicebox_v and bin/dicebox_1.0.1
#
# 2. make v            : Build DICEBOX-V OpenMP binary only
#                        -> Produces: bin/dicebox_v
#
# 3. make origin       : Build original DICEBOX-1.0.1 reference binary
#                        -> Produces: bin/dicebox_1.0.1
#
# 4. make clean        : Remove binary directory bin/ and temporary build files (*.o, *.mod)
#
# 5. make help         : Display this help message in the terminal
# ==============================================================================

FC          = gfortran
FFLAGS      = -fopenmp -O2 -fwrapv -std=legacy -w -ffp-contract=off
FFLAGS_ORIG = -O2 -std=legacy -w -ffp-contract=off

SRC_DIR     = src
ORIG_DIR    = v1.0.1
BIN_DIR     = bin

SRC_V       = $(SRC_DIR)/dicebox_v.for
SRC_ORIG    = $(ORIG_DIR)/dicebox_1.0.1.for

BIN_V       = $(BIN_DIR)/dicebox_v
BIN_ORIG    = $(BIN_DIR)/dicebox_1.0.1

.PHONY: all v origin clean help

all: $(BIN_V) $(BIN_ORIG)

v: $(BIN_V)

origin: $(BIN_ORIG)

$(BIN_V): $(SRC_V)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) $(SRC_V) -o $(BIN_V)
	@echo "[OK] Built $(BIN_V)"

$(BIN_ORIG): $(SRC_ORIG)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS_ORIG) $(SRC_ORIG) -o $(BIN_ORIG)
	@echo "[OK] Built $(BIN_ORIG)"

clean:
	rm -rf $(BIN_DIR)
	rm -f *.o *.mod
	@echo "[OK] Cleaned build artifacts."

help:
	@echo "=============================================================================="
	@echo "                      DICEBOX-V MAKEFILE USAGE GUIDE                          "
	@echo "=============================================================================="
	@echo "  make (or make all)  : Build both DICEBOX-V and DICEBOX-1.0.1 reference"
	@echo "  make v              : Build DICEBOX-V OpenMP only (bin/dicebox_v)"
	@echo "  make origin         : Build DICEBOX-1.0.1 reference only (bin/dicebox_1.0.1)"
	@echo "  make clean          : Clean bin/ and temporary build files"
	@echo "  make help           : Display this help message"
	@echo "=============================================================================="

