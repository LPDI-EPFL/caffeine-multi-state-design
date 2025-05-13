# caffeine-multi-state-design

**caffeine-multi-state-design** (varCOMETS) is a Python-based program developed by the Laboratory of Protein Design and Immunoengineering (LPDI) at EPFL for multi-state design. It focuses on multi-state protein design, aiming to optimize protein sequences that can adopt multiple conformations or functional states. This implementation was used to design **caffeine-inducible nanobody heterodimers** with minimized off-target homodimerization, enabling precise control of cellular signaling in synthetic biology applications.

If you use this code please cite our paper, and the original COMETS multi-state design paper [2].

![Caffeine multi state design](images/Caffeine.png)

## Features

- Implements multi-state design algorithms to identify sequences compatible with multiple protein states.
- Utilizes Python/pymol for scripting.
- Incorporates Rosetta tools for energy calculations, rotamer sampling, and structural analysis.

## Installation

To use this repository, clone it to your local machine:

```bash
git clone https://github.com/LPDI-EPFL/caffeine-multi-state-design.git
cd caffeine-multi-state-design
```

## Usage
The main script for running the multi-state design is design.py. To execute the design process:

```bash
python msd_caf.py
```

This will process the input structures and generate designed sequences compatible with the specified multiple states.

## Repository Structure

``design.py``: Main script to initiate the multi-state design process.

``msd.py, msd_caf.py``: Modules containing core functions for multi-state design algorithms.

``show_in_pymol.py``: Utility script to visualize structures and designs in PyMOL.

``test-mplp-orig.py, test_pyrosetta.py``: Test scripts for validating design methods.

``input/``: Directory containing input PDB files and configuration settings.

``output/``: Directory where output files, including designed sequences and structures, are saved.
GitHub

## License
This project is licensed under the MIT License. See the LICENSE file for details.

## References
If you use this code, please cite:

[1] Scheller L. et al. _"Humanized Caffeine-Inducible Systems for Controlling Cellular Functions"_, 2025

[2] Hallen M. & Donald B.R., _"COMETS (Constrained Optimization of Multistate Energies by Tree Search): A provable and efficient protein design algorithm to optimize binding affinity and specificity with respect to sequence."_ Journal of Computational Biology 23.5 (2016): 311-321.



