# mdspass2

`mdspass2` is an in-house molecular dynamics (MD) visualization and analysis tool developed in Y. Umeno Laboratory (Institute of Industrial Science, The University of Tokyo).  
It supports several empirical interatomic potentials and provides an interactive GUI for structure editing, MD simulation, visualization, NEB analysis, phonon analysis, and more.

This repository contains the source code, auxiliary libraries, sample files, and pre-built binaries for some platforms.

---

## Features

- Interactive MD visualization (GUI)
- Support for various empirical potentials  
  (Morse, EAM, ADP, Tersoff, AIREBO (torsional term excluded), etc.)
- NEB (Nudged Elastic Band) analysis
- Phonon dispersion calculation from small unit cells
- Structure editing (add/delete/move atoms, change species)
- Stress analysis (global & local)
- Periodic boundary condition handling
- Configuration measurement tools (distances, bond angles)
- Energy/force-only update mode (no MD)

---

## Directory Structure

```

mdspass2/            # Source code (C++) for Linux
├── CONFIG/          # Example configuration files
├── SETDAT/          # Example setting files
├── TUTORIAL/        # Tutorials
├── forMac/          # Source code, auxiliary libraries & build scripts for macOS
├── forWin/          # Pre-compiled binary for Windows
├── old/             # Old source code
├── pot/             # Potential files
├── sample_files/    # Sample CONFIG and SETDAT files
└── README.md        # This file

```

---

## Download and build
**NB: See License below before downloading.**

### Windows (pre-compiled binary)
- https://github.com/yoshiumeno/mdspass2/tree/main/forWin

### Linux/macOS (compile from source)
- Source code: 
 https://github.com/yoshiumeno/mdspass2 (for Linux)
 https://github.com/yoshiumeno/mdspass2/tree/main/forMac (for macOS)

### How to compile (Linux/Mac)

The detailed instruction is available on the laboratory website: https://www.cmsm.iis.u-tokyo.ac.jp/en/code_info_mdspass2_en.html (English)
 https://www.cmsm.iis.u-tokyo.ac.jp/en/code_info_mdspass2.html (Japanese)


---

## Known Issues

Please refer to the “Known problems” section on the laboratory website.
Typical examples include:

* Occasional segmentation faults on certain Linux graphics drivers
* GLUI-related GUI limitations
* macOS build issues depending on gcc/clang versions

---

## Update History

The detailed changelog is maintained on the laboratory website: https://www.cmsm.iis.u-tokyo.ac.jp/en/code_info_mdspass2_en.html (English)
https://www.cmsm.iis.u-tokyo.ac.jp/code_info_mdspass2.html (Japanese)

Major updates include:

* ReaxFF compatible version merge (2021)
* NEB improvements (2021)
* Bug fixes for Tersoff / ADP / EAM Mishin / GEAM (multiple versions)
* AIREBO implementation (2013)
* macOS support added (2013)

---

## License

At this moment, the use of MDSPASS2 is limited to our collaborative research projects. If you wish to download or use the software, please make sure to contact Umeno Laboratory in advance.

---

## Contact
Prof. Yoshitaka Umeno
Umeno Laboratory
Institute of Industrial Science, University of Tokyo

* Website: https://www.cmsm.iis.u-tokyo.ac.jp/index.html
* Email: umeno@iis.u-tokyo.ac.jp 


---
