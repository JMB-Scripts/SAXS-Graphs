# SAXS-Graphs
Regain control of your SAXS data.
This Python script is designed to process SAXS (Small Angle X-ray Scattering) data and to generate results such as the Guinier approximation, the Kratky graph and the volume of correlation in a text file format.
Then you can make YOUR OWN representations for these graphs, with YOUR preferred software (Excel, Prism, Origin etc.). 
## User manual

# Prerequisites
- Python 3.x
- Necessary packages: numpy, matplotlib, scipy
- SAXS Data need to be in Å-1

## Command syntax

python Sexier-vXX.py filename.dat qmin_offset qmax_offset

- `filename.dat` : the name of the .dat file containing the SAXS experimental data.
- `qmin_offset` : the offset (in a number of lines) to be added to the first usable line to determine qmin, use the value from PRIMUS (Range) or RAW (nmin).
- `qmax_offset`: the offset (in a number of lines) to be added to the first usable line to determine qmax use the value from PRIMUS (Range) or RAW (nmax).

## Features

 1. Guinier approximation :
 - Read .dat file and determine a first usable line.
 - Data extraction for q and I(q) in the selected range.
 - Linear regression to calculate Rg (radius of gyration) and I0 (intensity at q=0).
 - Write data to text file.
 - Display graph with experimental points and theoretical curve.

 2. Kratky 2:
 - Extract data in the selected range for Kratky.
 - Calculation and normalisation of values for Kratky (𝑞𝑅𝑔)^2.𝐼(𝑞)/𝐼(0) vs 𝑞𝑅𝑔.
 - Write data to a text file.
 - Display Kratky graph.

 3. Correlation volume (CV):
 - Extract data up to q=0.3 or up to 8/Rg.
 Calculate the integral of the product I(q)*q.
 - Calculation of VC, QR (quality factor) and MW (molecular weight).
 - Write data to text file.

 4. Summary file:
 - Write values for Rg, I0, qmin_Rg, qmax_Rg, MW to a text file.

## Output
 The script will generate the following files:
 - filename_rg.txt: Guinier approximation data (q^2, ln(I_exp), ln(I_theo), normalized residuals)
 - file_name_Normalized-Kratky.txt: data for normalized Kratky (x, y)
 - file_name_VC.txt: data for VC (q, I(q)*q)
 - file_name_Summary.txt: summary file (Rg, I0, qmin_Rg, qmax_Rg, MW)

![image](https://github.com/user-attachments/assets/e1155d9f-36c8-4d87-968b-0bbe50d49973)
   
## Plots:
  1- Guinier Fit (Rg, I(0), qmin*Rg qmax*Rg,nbeg, nend) 
  2- Residuals of the fit (check aggregation and or repulsion, here looks nice)
  3- Normalised Kratky plot (presence of disordered regions)
  4- Plot of the Volume of correlation and cumulative integral (evaluate the q max for MW determination)
  
![image](https://github.com/user-attachments/assets/a5c137c3-fd36-4dc0-9e65-fe48037712da)

## Notes:
1. It's possible to make an exe file for Windows using "pyinstaller", to distribute the script on computers without Python.
2. I can also provide the stand-alone version for Windows upon request. 
