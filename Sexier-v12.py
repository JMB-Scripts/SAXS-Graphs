# User manual
# Jean-Marie Bourhis (jean-marie.bourhis@ibs.fr) and GPT via mac-gpt (because I'm a "tanche" (kind of fish) with python)
# 2023
# This Python script is designed to asses the quality of your SAXS data :
# 1- From Guinier approximation, I(0), Rg
# 2- the Kratky graph (presence or not of disordered region)
# 3- the correlation volume to get the MW.
#
# It will allow you to keep all these information in several simple text files, allowing you to use your own graphical software (Excel, Graphpad, Origin, ...).
# 
## Prerequisites
#- Python 3.x
#- Necessary packages: numpy, matplotlib, scipy
#- SAXS Data need to be in Å-1

## Command syntax
# python script.py filename.dat qmin_offset qmax_offset
#- `filename.dat` : the name of the .dat file containing the SAXS experimental data.
#- `qmin_offset` : the offset (in number of lines) to be added to the first usable line to determine qmin, use the value from PRIMUS or RAW.
#- `qmax_offset`: the offset (in number of lines) to be added to the first usable line to determine qmax use the value from PRIMUS or RAW.
## Features
# 1. Guinier approximation :
# - Read .dat file and determine first usable line.
# - Extraction of data for q and I(q) in the selected range.
# - Linear regression to calculate Rg (radius of gyration) and I0 (intensity at q=0).
# - Write data to text file.
# - Display graph with experimental points and theoretical curve.
# 2. Kratky 2:
# - Extract data in selected range for Kratky.
# - Calculation and normalization of values for Kratky (𝑞𝑅𝑔)^2.𝐼(𝑞)/𝐼(0) vs 𝑞𝑅𝑔.
# - Write data to text file.
# - Display Kratky graph.
# 3. Correlation volume (CV):
# - Extract data up to q=0.3.
# Calculate the integral of the product I(q)*q.
# - Calculation of VC, QR (quality factor) and MW (molecular weight).
# - Write data to text file.
# 4. Summary file:
# - Write values for Rg, I0, qmin_Rg, qmax_Rg, MW to a text file.
## Output
# The script will generate the following files:
#   - filename_rg.txt: Guinier approximation data (q^2, ln(I_exp), ln(I_theo), normalised residuals)
#   - file_name_Normalized-Kratky.txt: data for normalised Kratky (x, y)
#   - nom_du_fichier_VC.txt: data for VC (q, I(q)*q)
#   - file_name_Summary.txt: summary file (Rg, I0, qmin_Rg, qmax_Rg, MW)
 
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# Check that arguments are supplied correctly
if len(sys.argv) < 4:
    print("Veuillez fournir les arguments nécessaires : nom_du_script.py nom_du_fichier.dat qmin_offset qmax_offset")
    sys.exit(1)

data_file = str(sys.argv[1])
qmin_offset = int(sys.argv[2])
qmax_offset = int(sys.argv[3])

# Extract the first 5 characters of the source file name
source_file_prefix = os.path.splitext(os.path.basename(sys.argv[1]))[0]

# Define output file names with prefixes
output_file = f'{source_file_prefix}_01_Rg.txt'
output_file2 = f'{source_file_prefix}_02_Norm-Kratky.txt'
output_file3 = f'{source_file_prefix}_03_VC.txt'
output_file4 = f'{source_file_prefix}_04_Summary.txt'
#output_file5 = f'5-{source_file_prefix}_Integrale.txt'

#####################################
#
# Reading the file
#
#####################################

# Read the file and find the first usable line

first_usable_line = None
with open(data_file, 'r') as file:
    for line_number, line in enumerate(file):
        line = line.strip()
        if line and len(line.split()) >= 3:
            try:
                values = [float(x) for x in line.split()[:3]]
                if not any(np.isnan(values)):
                    first_usable_line = line_number
                    break
            except ValueError:
                continue

        # If the first usable line has already been found, exit the loop.
        if first_usable_line is not None:
            break

if first_usable_line is None:
    print("Aucune ligne utilisable trouvée dans le fichier.")
else:
    print("Première ligne utilisable :", first_usable_line)

# Lecture du fichier et recherche de la dernière ligne utilisable
last_usable_line = None
with open(data_file, 'r') as file:
    lines = file.readlines()
    for line_number, line in reversed(list(enumerate(lines, start=1))):
        line = line.strip()
        if line and len(line.split()) >= 3:
            try:
                values = [float(x) for x in line.split()[:3]]
                if not any(np.isnan(values)):
                    last_usable_line = line_number
                    break
            except ValueError:
                continue
if last_usable_line is None:
    print("Aucune ligne utilisable trouvée dans le fichier.")
else:
    print("Numéro de la dernière ligne utilisable :", last_usable_line)

nbr_lines = last_usable_line - first_usable_line 
print("nbr of lines", nbr_lines)

#####################################
#
# Guinier
#
#####################################
    # Load data from .dat file

data = np.loadtxt(data_file, skiprows=first_usable_line, max_rows=nbr_lines)


   # Extract data columns for q and I(q)
q = data[:, 0]

    # Find the index of the column containing the values of I
I = data[:, 1]
    # Find the index of the column containing the errors on I
E = data[:, 2]
    # Select q values corresponding to qmin and qmax
qmin = q[qmin_offset]
qmax = q[qmax_offset]
I_qmin = I[qmin_offset]
I_qmax = I[qmax_offset]
    
    # Select q range for linear regression
q_range = q[(q >= qmin) & (q <= qmax)]
I_range = I[(q >= qmin) & (q <= qmax)]
E_range = E[(q >= qmin) & (q <= qmax)]
 # Perform linear regression to find Rg and I0
m, c = np.polyfit(q_range**2, np.log(I_range), deg=1)

Rg = np.sqrt(-3 * m)
I0 = np.exp(c)

 # Calculate products qmin * Rg and qmax * Rg
qmin_Rg = qmin * Rg
qmax_Rg = qmax * Rg

 # Display results
print("Rg:", Rg)
print("I0:", I0)
print("qmin * Rg:", qmin_Rg)
print("qmax * Rg:", qmax_Rg)

   # Calculate q^2
q_carre =  q_range * q_range
   #  Calculate theoretical values based on linear regression
ln_intensity_theoretical = m * q_range**2 + c
   # Calculate ln of I_range 
ln_intensity_exp = np.log(I_range)
   # Calculate normalized residuals
residuals = (np.log(I_range) - ln_intensity_theoretical) / E_range

# Create a text file containing the data
with open(output_file, 'w') as file:
        file.write("q**2\tln(I_exp)\tln(I_theo)\tnormalized-residuals\n")
        for i in range(len(q_range)):
            file.write(f"{q_carre[i]}\t{ln_intensity_exp[i]}\t{ln_intensity_theoretical[i]}\t{residuals[i]}\n")
    #print("Données enregistrées dans", output_file)

#####################################
#
#  Kratky
#
#####################################
#
# Creation of a txt file for a Kratky normalized by the Rg
# (𝑞𝑅𝑔)^2.𝐼(𝑞)/𝐼(0) vs 𝑞𝑅𝑔
#
###################

#

firstqmin = q[1]
firstqmax = q[nbr_lines - 1]

#  range definition
q_full =  q[(q >= firstqmin ) & (q <= firstqmax)]
I_full =  I[(q >= firstqmin ) & (q <= firstqmax)]

# Calculate qrg
xkrat = q_full * Rg
# Calculate I(q)/I(0) named nI
zkrat = xkrat * xkrat * I_full
ykrat = np.divide(zkrat,I0)
# Create a text file containing data
with open(output_file2, 'w') as file:
        file.write("q.Rg\t(q.Rg)^2.(I(q)/I(0)\n")
        for j in range(len(q_full)):
            file.write(f"{xkrat[j]}\t{ykrat[j]}\n")



#####################################
#
# Volume of correlation 
#
#####################################

# Filter data up to q=0.3
q_filtered = q[q <= 0.3]
I_filtered = I[q <= 0.3]

# define the Y axis with yvc
yvc = I_filtered * q_filtered

# define the integrale Intgr
Intgr = integrate.simps(yvc, q_filtered)

# calculate VC, QR and MW
VC = I0 / Intgr
QR = VC**2 / Rg
MW1 = QR / 0.1231
print("MW cut q=0.3:", MW1)

# Filter data up to q=8/rg
q_vc = 8 / Rg
q_alt = q_full[q_full <= q_vc]
I_alt = I_full[q_full <= q_vc]

# Define Y for vc yvc
yvc_alt = I_alt * q_alt
# Define the integrale named Intgr_alt
Intgr_alt = integrate.simps(yvc_alt, q_alt)

# Calculate VC, QR and MW
VC_alt = I0 / Intgr_alt
QR_alt = VC_alt**2 / Rg
MW_alt = QR_alt / 0.1231
print("MW cut q=8/Rg:", MW_alt)

# For the plot

# Define Yvc_plot
yvc_plot = I_full * q_full
# Define the integrale named Intgr_alt
Intgr_plot = integrate.simps(yvc_plot, q_full)
# Output file with I(q)q vs q et Intgr vs q

with open(output_file3, 'w') as file:
        file.write("q\tI(q)*q\n")
        for i in range(len(q_full)):
            file.write(f"{q_full[i]}\t{yvc_plot[i]}\n")

#####################################
#
# Summary
#
#####################################
# Output summary 
with open(output_file4, 'w') as file:
    file.write("Rg\n")
    file.write(f"{Rg}\n")
    file.write("I0\n")
    file.write(f"{I0}\n")
    file.write("qmin_Rg\n")
    file.write(f"{qmin_Rg}\n")
    file.write("qmax_Rg\n")
    file.write(f"{qmax_Rg}\n")
    file.write("MW from Vc cut q=0.3\n")
    file.write(f"{MW1}\n")
    file.write("MW from Vc cut q=8/Rg\n")
    file.write(f"{MW_alt}\n")

#  Display Guinier graph with experimental points and theoretical curve
plt.figure(figsize=(10, 6))
plt.errorbar(q_range**2, np.log(I_range), yerr=E_range, fmt='o', markersize=4, label='expériment')
plt.plot(q_range**2, ln_intensity_theoretical, label='Fit')
plt.xlabel('q**2')
plt.ylabel('ln(I(q))')
plt.legend()
plt.title('Guinier')
plt.grid(True)
plt.show()

# Display graph of standardized residuals
plt.figure(figsize=(10, 4))
plt.plot(q_range, residuals, 'o', markersize=4)
plt.xlabel('q**2')
plt.ylabel('Guinier normalised residuals')
plt.title('Guinier normalised residuals')
plt.grid(True)
plt.show()

#  Draw the graph Kratky normalised
plt.plot(xkrat, ykrat)
plt.xlabel('(qRg)^2')
plt.ylabel('(qRg)^2 * I(q) / I(0)')
plt.title('Normalised Kratky Plot')
plt.grid(True)
plt.show()

#  Draw the graph VC
plt.plot(q_full, yvc_plot)
plt.xlabel('(qRg)^2')
plt.ylabel('(qRg)^2 * I(q) / I(0)')
plt.title('Volume of Correlation')
plt.grid(True)
plt.show()
