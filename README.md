# NOMAD Molecular Parser & Normalizer
# nomad-plugin-molecules

This plugin provides parsing and normalization for molecular structure data within the NOMAD framework.

## Installation
```
pip install .
```

### Installing OpenBabel via pip
OpenBabel is a key dependency for this project and is included in `pyproject.toml` and `requirements.txt` as `openbabel-wheel`, so it should be installed automatically. If for any reason it isn’t, you can manually install it with:
```sh
pip install openbabel-wheel
```
**Note:**
OpenBabel relies on system libraries such as libxrender1 and libxext6. In minimal installations (including Docker containers) these libraries may be missing. If you encounter errors indicating that `libXrender.so.1` or `libXext.so.6` is missing, you’ll need to install these libraries manually.

After installing OpenBabel, verify the installation by running:

```sh
obabel -V
```
This command should display the version of Open Babel, confirming that it is installed correctly.

#### Installing System Dependencies

**Debian/Ubuntu:**
```sh
sudo apt-get update
sudo apt-get install libxrender1 libxext6
```

**Fedora/CentOS/RHEL:**
```sh
sudo dnf install libXrender libXext
```

**Arch Linux:**
```sh
sudo pacman -S libxrender libxext
```

**macOS:**
Install XQuartz to provide the necessary X11 libraries.

**Windows:**
If you are running a Linux environment (e.g., via WSL or Cygwin), use the appropriate Linux commands. Otherwise, Windows typically does not require these libraries unless you are running Linux-based tools.

## Usage
```

```

## Features
- Converts structures to InChIKeys.
- Queries PubChem database for additional molecular details.