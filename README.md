# caffeine-multi-state-design

**caffeine-multi-state-design** is a Python-based framework developed by the Laboratory of Protein Design and Immunoengineering (LPDI) at EPFL for multi-state desing. It focuses on multi-state protein design, aiming to optimize protein sequences that can adopt multiple conformations or functional states.

## Features

- Implements multi-state design algorithms to identify sequences compatible with multiple protein states.
- Utilizes Python for scripting and automation of design workflows.
- Incorporates tools for energy calculations, rotamer sampling, and structural analysis.

## Installation

To use this repository, clone it to your local machine:

```bash
git clone https://github.com/LPDI-EPFL/caffeine-multi-state-design.git
cd caffeine-multi-state-design```
Installation
To use this repository, clone it to your local machine:

bash
Copy
Edit
git clone https://github.com/LPDI-EPFL/caffeine-multi-state-design.git
cd caffeine-multi-state-design
Install the required Python packages:

bash
Copy
Edit
pip install -r requirements.txt
Note: Ensure you have Python 3.7 or higher installed.
GitHub

## Usage
The main script for running the multi-state design is design.py. To execute the design process:

```
bash
python design.py
```

This will process the input structures and generate designed sequences compatible with the specified multiple states.


Repository Structure

``design.py``: Main script to initiate the multi-state design process.

msd.py, msd_caf.py: Modules containing core functions for multi-state design algorithms.

show_in_pymol.py: Utility script to visualize structures and designs in PyMOL.

test-mplp-orig.py, test_pyrosetta.py: Test scripts for validating design methods.

input/: Directory containing input PDB files and configuration settings.

output/: Directory where output files, including designed sequences and structures, are saved.
GitHub

## License
This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgments
This tool was developed by the Laboratory of Protein Design and Immunoengineering (LPDI) at École Polytechnique Fédérale de Lausanne (EPFL). For more information, visit the LPDI GitHub page.

