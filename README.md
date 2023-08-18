# SAXS-Graphs
Regain control of your SAXS data
 Jean-Marie Bourhis (jean-marie.bourhis@ibs.fr) and GPT via mac-gpt (because I'm a "tanche" (kind of fish) with python)
 2023
 This Python script is designed to process SAXS (Small Angle X-ray Scattering) data and to generate results such as the Guinier approximation, the Kratky graph and the correlation volume.

## User manual

# Prerequisites
- Python 3.x
- Necessary packages: numpy, matplotlib, scipy
- SAXS Data need to be in Å-1

## Command syntax

python Sexier-vXX.py filename.dat qmin_offset qmax_offset

- `filename.dat` : the name of the .dat file containing the SAXS experimental data.
- `qmin_offset` : the offset (in number of lines) to be added to the first usable line to determine qmin, use the value from PRIMUS or RAW.
- `qmax_offset`: the offset (in number of lines) to be added to the first usable line to determine qmax use the value from PRIMUS or RAW.

## Features

 1. Guinier approximation :
 - Read .dat file and determine first usable line.
 - Extraction of data for q and I(q) in the selected range.
 - Linear regression to calculate Rg (radius of gyration) and I0 (intensity at q=0).
 - Write data to text file.
 - Display graph with experimental points and theoretical curve.

 2. Kratky 2:
 - Extract data in selected range for Kratky.
 - Calculation and normalization of values for Kratky (𝑞𝑅𝑔)^2.𝐼(𝑞)/𝐼(0) vs 𝑞𝑅𝑔.
 - Write data to text file.
 - Display Kratky graph.

 3. Correlation volume (CV):
 - Extract data up to q=0.3.
 Calculate the integral of the product I(q)*q.
 - Calculation of VC, QR (quality factor) and MW (molecular weight).
 - Write data to text file.

 4. Summary file:
 - Write values for Rg, I0, qmin_Rg, qmax_Rg, MW to a text file.

## Output
 The script will generate the following files:
 - filename_rg.txt: Guinier approximation data (q^2, ln(I_exp), ln(I_theo), normalized residuals)
 - file_name_Normalized-Kratky.txt: data for normalized Kratky (x, y)
 - nom_du_fichier_VC.txt: data for VC (q, I(q)*q)
 - file_name_Summary.txt: summary file (Rg, I0, qmin_Rg, qmax_Rg, MW)
