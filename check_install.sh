#!/usr/bin/env bash
# ==============================================================================
# check_install.sh — Installation Verification Script for DICEBOX-V
# Computer Physics Communications Program Library
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${SCRIPT_DIR}/bin/dicebox_v"
SRC="${SCRIPT_DIR}/src/dicebox_v.for"
EX1="${SCRIPT_DIR}/examples/Example1"

echo "================================================================================"
echo "    DICEBOX-V: Installation Verification Gate"
echo "================================================================================"

# Check for compiler
if ! command -v gfortran &> /dev/null; then
    echo "[ERROR] gfortran compiler not found in PATH."
    echo "        Please install GNU Fortran with OpenMP 5.0 support (verified with gfortran >= 13)."
    exit 1
fi

# Build binary if not already present
if [ ! -f "${BIN}" ]; then
    echo "[INFO] Compiling DICEBOX-V binary via Makefile..."
    if command -v make &> /dev/null; then
        make -C "${SCRIPT_DIR}" v
    else
        mkdir -p "${SCRIPT_DIR}/bin"
        gfortran -fopenmp -O2 -fwrapv -std=legacy -w -ffp-contract=off \
            "${SRC}" -o "${BIN}"
    fi
fi

if [ ! -x "${BIN}" ]; then
    echo "[ERROR] Binary ${BIN} could not be built or is not executable."
    exit 1
fi

echo "[INFO] Binary verified: ${BIN}"

# Create isolated temporary directory for test runs
TEST_DIR=$(mktemp -d /tmp/dicebox_test_install.XXXXXX)
trap 'rm -rf "${TEST_DIR}"' EXIT INT TERM

echo "[INFO] Running installation test in temporary sandbox: ${TEST_DIR}"
mkdir -p "${TEST_DIR}/run_test"
find "${EX1}" -maxdepth 1 -type f -exec cp {} "${TEST_DIR}/run_test/" \;
rm -f "${TEST_DIR}/run_test"/DICE.PRO "${TEST_DIR}/run_test"/PSF_GS.DAT "${TEST_DIR}/run_test"/PSF_INI.DAT "${TEST_DIR}/run_test"/EVENTS.* "${TEST_DIR}/run_test"/LEVELS.*

NTH_CHECK="${OMP_NUM_THREADS:-2}"
if [ "${NTH_CHECK}" -lt 1 ]; then
    NTH_CHECK=2
fi
echo "[INFO] Executing Example 1 with OMP_NUM_THREADS=${NTH_CHECK}..."
(
    cd "${TEST_DIR}/run_test"
    OMP_NUM_THREADS="${NTH_CHECK}" "${BIN}" > run.log 2>&1
) || {
    EXIT_CODE=$?
    echo "[ERROR] DICEBOX-V (OMP_NUM_THREADS=${NTH_CHECK}) exited with non-zero status: ${EXIT_CODE}"
    echo "--- Last 20 lines of log: ---"
    tail -n 20 "${TEST_DIR}/run_test/run.log"
    exit 1
}

MISSING=0
EXPECTED_FILES=(
    DICE.PRO PSF_GS.DAT PSF_INI.DAT
    LEVELS.001 LEVELS.002 LEVELS.003
    EVENTS.S001.R001 EVENTS.S001.R002
    EVENTS.S002.R001 EVENTS.S002.R002
    EVENTS.S003.R001 EVENTS.S003.R002
)
for expected_file in "${EXPECTED_FILES[@]}"; do
    if [ ! -s "${TEST_DIR}/run_test/${expected_file}" ]; then
        echo "[ERROR] Expected output file '${expected_file}' was not generated or is empty."
        MISSING=1
    fi
done

if [ ${MISSING} -ne 0 ]; then
    echo "[FAIL] Installation verification failed."
    exit 1
fi

echo "[INFO] All ${#EXPECTED_FILES[@]} expected Example1 output files were generated successfully."

echo "================================================================================"
echo "    [PASS] DICEBOX-V installed and verified successfully!"
echo "================================================================================"
exit 0
