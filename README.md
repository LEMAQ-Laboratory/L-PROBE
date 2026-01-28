# L-PROBE: Local PROminence-Based Elimination

**Local PROminence-Based Elimination (L-PROBE): A fast and robust peak selection algorithm for Raman spectroscopy data reduction and chemometric analysis.**

## Overview
L-PROBE is a chemometric algorithm designed to optimize variable selection in spectroscopic data, specifically focused on Raman spectra. By utilizing local prominence criteria and neighborhood filtering, it identifies chemically relevant variables, significantly reducing data dimensionality while preserving essential information for multivariate modeling.

## Key Features
* **Efficient Data Reduction:** Drastically reduces the number of variables for faster computation.
* **Robust Peak Selection:** Focuses on signals with high local prominence, minimizing noise interference.
* **Automated Workflow:** Reduces the need for manual feature selection in large datasets.

## Repository Structure
* `lprobe.m`: Core MATLAB function implementing the L-PROBE algorithm.
* `.gitignore`: Standard filters for MATLAB temporary files.

## Installation
To use L-PROBE, simply download the `lprobe.m` file and add it to your MATLAB path, or clone this repository:
```bash
git clone https://github.com/LEMAQ-Laboratory/L-PROBE.git
