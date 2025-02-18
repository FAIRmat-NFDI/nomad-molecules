#!/bin/bash

set -e  # Exit on error
SCRIPT_DIR=$(dirname "$(realpath "$0")")  # Get the directory of the script

echo "Updating package list..."
apt update

echo "Installing required system dependencies..."
apt install -y cmake g++ libeigen3-dev libboost-all-dev libxml2-dev libz-dev swig python3-dev python3-pip git

echo "Installing Open Babel..."
git clone https://github.com/openbabel/openbabel.git
cd openbabel
mkdir build && cd build
cmake .. -DPYTHON_BINDINGS=ON -DPYTHON_EXECUTABLE=$(which python3) -DRUN_SWIG=ON
make -j$(nproc)
make install  # Ensure system-wide installation

cd "$SCRIPT_DIR"

echo "Setting up Open Babel environment..."
echo "/usr/local/lib" | tee -a /etc/ld.so.conf.d/openbabel.conf
ldconfig
echo 'export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc
source ~/.bashrc
ldconfig

echo "Verifying Open Babel installation..."
obabel -H || { echo "Open Babel installation failed"; exit 1; }

echo "Installing molid..."
# Ensure pip is up-to-date
pip install --upgrade pip setuptools wheel
# Install molid from its repository
pip install git+https://gitlab.mpcdf.mpg.de/tdenell/molecule-identification.git

echo "Verifying molid installation..."
python3 -c "import molid; print(molid.__version__)" || { echo "molid installation failed"; exit 1; }

echo "Installing your normalizer (nomad-plugin-molecules)..."
# Navigate to the directory containing pyproject.toml
cd "$SCRIPT_DIR"
pip install .

echo "Verifying normalizer installation..."
python3 -c "import nomad_plugin_molecules; print('Normalizer installed successfully')" || { echo "Normalizer installation failed"; exit 1; }

echo "Installation complete!"
