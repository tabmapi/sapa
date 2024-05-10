# SAPA
This repository contains the files necessary for running symmetry-adapted PDF analysis with Topas v6 and v7.

## Requirements
* Topas v6 or v7
* Python 3.6 or later
* HDF5
* Topas_for_PDF macros - available at github.com/pachater/topas

## Installation

Currently, SAPA is not available as a package. To use, download sapa.py and sapa_utils_hdf.py into a folder on your computer. Add this folder to your python path, which you can do by searching 'Edit System Environment Variables' in the taskbar, clicking the Environment Variables button, and under User variables, add a new variable with name PYTHONPATH, with the value as the folder you placed the downloaded files. Your python distribution will now know to look in that folder to import files (in addition to default locations), so you can import sapa from anywhere.

## Authors

This code is predominantly written and maintained by Tobias Bird, with contributions from Mark Senn.

## Contact

For any issues with the code, contact tobias (dot) bird (at) diamond.ac.uk


