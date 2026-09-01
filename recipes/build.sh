#!/bin/bash
set -exuo pipefail

SHARE_DIR="${PREFIX}/share/spechla"

# --- Build SpecHap ---
mkdir -p bin/SpecHap/build
cd bin/SpecHap/build
cmake .. \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DCMAKE_BUILD_TYPE=Release
make -j${CPU_COUNT}
make install
cd "${SRC_DIR}"

# --- Build ExtractHAIRs ---
mkdir -p bin/extractHairs/build
cd bin/extractHairs/build
cmake .. \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DCMAKE_BUILD_TYPE=Release
make -j${CPU_COUNT}
make install
cd "${SRC_DIR}"

# --- Install scripts to share dir ---
mkdir -p "${SHARE_DIR}/script"
cp -r script/* "${SHARE_DIR}/script/"

# --- Install reference construction assets ---
mkdir -p "${SHARE_DIR}/reference_assets"
cp -r share/reference_assets/* "${SHARE_DIR}/reference_assets/"

# --- Install utility scripts to bin ---
install -m 755 bin/blast2sam.pl "${PREFIX}/bin/"
install -m 755 bin/vcf-combine.py "${PREFIX}/bin/"

# --- Install wrapper: spechla ---
cat > "${PREFIX}/bin/spechla" << 'EOF'
#!/bin/bash
set -euo pipefail
export SPECHLA_DB="${SPECHLA_DB:-${CONDA_PREFIX}/share/spechla/db}"
export SPECHLA_SCRIPT="${SPECHLA_SCRIPT:-${CONDA_PREFIX}/share/spechla/script}"
exec bash "${SPECHLA_SCRIPT}/whole/SpecHLA.sh" "$@"
EOF
chmod +x "${PREFIX}/bin/spechla"

# --- Install wrapper: spechla-extract-hla-reads ---
cat > "${PREFIX}/bin/spechla-extract-hla-reads" << 'EOF'
#!/bin/bash
set -euo pipefail
export SPECHLA_DB="${SPECHLA_DB:-${CONDA_PREFIX}/share/spechla/db}"
export SPECHLA_SCRIPT="${SPECHLA_SCRIPT:-${CONDA_PREFIX}/share/spechla/script}"
exec bash "${SPECHLA_SCRIPT}/ExtractHLAread.sh" "$@"
EOF
chmod +x "${PREFIX}/bin/spechla-extract-hla-reads"

# --- Install wrapper: spechla-long-read ---
cat > "${PREFIX}/bin/spechla-long-read" << 'EOF'
#!/bin/bash
set -euo pipefail
export SPECHLA_DB="${SPECHLA_DB:-${CONDA_PREFIX}/share/spechla/db}"
export SPECHLA_SCRIPT="${SPECHLA_SCRIPT:-${CONDA_PREFIX}/share/spechla/script}"
exec python3 "${SPECHLA_SCRIPT}/long_read_typing.py" "$@"
EOF
chmod +x "${PREFIX}/bin/spechla-long-read"

# --- Install wrapper: spechla-assembly ---
cat > "${PREFIX}/bin/spechla-assembly" << 'EOF'
#!/bin/bash
set -euo pipefail
export SPECHLA_DB="${SPECHLA_DB:-${CONDA_PREFIX}/share/spechla/db}"
export SPECHLA_SCRIPT="${SPECHLA_SCRIPT:-${CONDA_PREFIX}/share/spechla/script}"
exec python3 "${SPECHLA_SCRIPT}/typing_from_assembly.py" "$@"
EOF
chmod +x "${PREFIX}/bin/spechla-assembly"

# --- Install wrapper: spechla-loh ---
cat > "${PREFIX}/bin/spechla-loh" << 'EOF'
#!/bin/bash
set -euo pipefail
export SPECHLA_DB="${SPECHLA_DB:-${CONDA_PREFIX}/share/spechla/db}"
export SPECHLA_SCRIPT="${SPECHLA_SCRIPT:-${CONDA_PREFIX}/share/spechla/script}"
exec perl "${SPECHLA_SCRIPT}/cal.hla.copy.pl" "$@"
EOF
chmod +x "${PREFIX}/bin/spechla-loh"

# --- Install wrapper: spechla-build-reference ---
cat > "${PREFIX}/bin/spechla-build-reference" << 'EOF'
#!/bin/bash
set -euo pipefail
export SPECHLA_SCRIPT="${SPECHLA_SCRIPT:-${CONDA_PREFIX}/share/spechla/script}"
export SPECHLA_ASSETS="${SPECHLA_ASSETS:-${CONDA_PREFIX}/share/spechla/reference_assets}"
exec python3 "${SPECHLA_SCRIPT}/build_reference.py" "$@"
EOF
chmod +x "${PREFIX}/bin/spechla-build-reference"
