# nomad-molecules

**NOMAD Molecular Normalizer**

A plugin for the [NOMAD](https://nomad-coe.eu/) framework that extracts, parses, and normalizes molecular structure data (e.g., small molecules) from computational archives, enriching them with cheminformatics metadata from an offline PubChem database.

---

## Features

* **Structure Extraction**: Identify and unwrap 0D molecular topologies from NOMAD archives.
* **InChIKey Conversion**: Compute InChIKeys for ASE `Atoms` objects using Open Babel.
* **Offline PubChem Lookup**: Query a local PubChem SQLite database for full or skeletal (InChIKey14) matches.
* **Metadata Injection**: Attach `smiles`, `inchi`, `inchi_key`, and `match_type` into `topology.cheminformatics`.
* **Configurable Limits**: Filter by minimum/maximum atom counts to skip very small or very large systems.
* **Extensible**: Easily override database paths and search modes via environment variables or Pydantic settings.

---

## Installation

### 1. Prerequisites

* Python ≥ 3.8
* ASE (Atomic Simulation Environment)
* NOMAD Python client (`nomad-client`)
* Open Babel

### 2. Install via pip

```sh
pip install nomad-plugin-molecules
```

> **Note:** The package depends on `openbabel-wheel` to provide Open Babel bindings.

### 3. System Libraries

Open Babel requires X11 rendering libraries on Linux. Install if missing:

* **Debian/Ubuntu:**

  ```sh
  sudo apt-get update && sudo apt-get install libxrender1 libxext6
  ```
* **Fedora/CentOS/RHEL:**

  ```sh
  sudo dnf install libXrender libXext
  ```
* **Arch Linux:**

  ```sh
  sudo pacman -S libxrender libxext
  ```
* **macOS:**

  * Install [XQuartz](https://www.xquartz.org/).
* **Windows:**

  * When using WSL or Cygwin, use the corresponding Linux commands above. Otherwise, Windows typically does not require these libraries unless you are running Linux-based tools.

Verify Open Babel:

```sh
obabel -V
```

---

## License

This project is licensed under the [BSD 3-Clause License](LICENSE).
