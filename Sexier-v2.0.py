
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
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
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
output_file4 = f'{source_file_prefix}_04_VC_integral.txt'
output_file5 = f'{source_file_prefix}_05_Summary.txt'
output_file6 = f'{source_file_prefix}_Graphs.png'
output_file7 = f'{source_file_prefix}_Graphs.svg'

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
#print("Rg:", Rg)
#print("I0:", I0)
#print("qmin * Rg:", qmin_Rg)
#print("qmax * Rg:", qmax_Rg)

   # Calculate q^2
q_carre =  q_range * q_range
   #  Calculate theoretical values based on linear regression
ln_intensity_theoretical = m * q_range**2 + c
   # Calculate ln of I_range 
ln_intensity_exp = np.log(I_range)
   # Calculate normalized residuals
residuals = (np.log(I_range) - ln_intensity_theoretical) / (E_range * 10)

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
#print("MW cut q=0.3:", MW1)

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
#print("MW cut q=8/Rg:", MW_alt)

##########
# Calculate cumulative integral of q * I(q) with respect to q
q_integral_cumulative = [integrate.simps(I_full[:i] * q_full[:i], q_full[:i]) for i in range(1, len(q_full) + 1)]

# Save cumulative integral values to a text file
with open(output_file4, 'w') as file:
    file.write("q\tCumulative Integral\n")  # Header
    for i in range(len(q_full)):
        file.write(f"{q_full[i]}\t{q_integral_cumulative[i]:.6f}\n")  # Save q and cumulative integral

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
# Volume of Porod using SAXSMOW approach
# V'=A+B*V
#
#####################################
# Filter data up to q=8/rg
q_porod = q_full[q_full <= q_vc]
I_porod = I_full[q_full <= q_vc]

# Calculate the integral Q_p of q^2 * I(q) with respect to q
q_squared_I = q_porod**2 * I_porod
Q_p = integrate.simps(q_squared_I, q_porod)
# Calculate Porod volume V_p
V_p = (2 * np.pi**2 * I0) / Q_p

# Calculate V'
# V'=A+B.V
A=(-2.114e6*(q_vc**4)+ 2.920e6*(q_vc**3)-1.472e6*(q_vc**2)+3.349e5*(q_vc)-3.577e4)
B=(12.09*(q_vc**3)-9.39*(q_vc**2)+3.03*(q_vc)+0.29)
V_prime= A+V_p*B
# Calculate MW 
MV_P =(V_prime/1.2)

#####################################
#
# Summary
#
#####################################
# Output summary 
with open(output_file5, 'w') as file:
    file.write("Rg\n")
    file.write(f"{Rg:.3f}\n")
    file.write("I0\n")
    file.write(f"{I0:.3f}\n")
    file.write("qmin_Rg\n")
    file.write(f"{qmin_Rg:.3f}\n")
    file.write("qmax_Rg\n")
    file.write(f"{qmax_Rg:.3f}\n")
    file.write("nbeg\n")
    file.write(f"{qmin_offset}\n")
    file.write("nend\n")
    file.write(f"{qmax_offset}\n")  
    file.write("MW from Vc cut q=0.3\n")
    file.write(f"{MW1:3f}\n")
    file.write("MW from Vc cut q=8/Rg\n")
    file.write(f"{MW_alt:3f}\n")
    file.write("Vp Porod\n")
    file.write(f"{V_p:2f}\n")
    file.write("MW Porod\n")
    file.write(f"{MV_P:2f}\n")
#####################################
#
# PLOT
#
#####################################

col1='#0D92F4'#plot guinier and residual
col2='#77CDFF'
col3='#F95454'
col4='#C62E2E'
col5='#F3F3E0'#boxbackground
col6='#BC7C7C' #cumulative 
# Create a figure with 2x2 subplots, four panels A,B,C,D
fig, axs = plt.subplots(2, 2, figsize=(13, 11))

# A- Upper-left: Form Factor plot (log(I) vs q)
# Filter out negative I_full values
positive_indices = I_full > 0
q_full_positive = q_full[positive_indices]
I_full_positive = I_full[positive_indices]
# Plot only positive values
axs[0, 0].plot(q_full_positive, np.log(I_full_positive), label='Log(I) vs q', color=col1)
axs[0, 0].set_xlabel('q')
axs[0, 0].set_ylabel('log(I(q))')
axs[0, 0].set_title(str(source_file_prefix) +' Form Factor' ) #add name of the file
axs[0, 0].legend()
axs[0, 0].grid(True)

# B-1-Lower-left: Guinier plot
axs[1, 0].errorbar(q_range**2, np.log(I_range), yerr=E_range, fmt='o', markersize=4, label='Experimental', color=col1)
axs[1, 0].plot(q_range**2, ln_intensity_theoretical, label='Fit', color=col4)
axs[1, 0].set_xlabel('q**2')
axs[1, 0].set_ylabel('ln(I(q))')
axs[1, 0].set_title('Guinier Plot')
axs[1, 0].ticklabel_format(axis="x", style="sci", scilimits=(0,0), useMathText=True) #sci test
axs[1, 0].legend()
axs[1, 0].grid(True)

# Add Guinier details as text
guinier_text = f'Rg = {Rg:.3f}\nI0 = {I0:.3f}\nqmin * Rg = {qmin_Rg:.3f}\nqmax * Rg = {qmax_Rg:.3f}\nnbeg = {qmin_offset}\nnend = {qmax_offset}'
axs[1, 0].text(0.72, 0.6, guinier_text, transform=axs[1, 0].transAxes, fontsize=10,
               verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor=col5, alpha=0.5))# Adjusted position and size: left, bottom,

# B-2-Create inset for residuals (smaller and repositioned)
inset_residuals = fig.add_axes([0.09, 0.1, 0.21, 0.09])  # Adjusted position and size: left, bottom, width, height  
inset_residuals.plot(q_range**2, residuals, 'o', markersize=3, color=col1, label='Residuals')
inset_residuals.axhline(y=0, color=col4, linewidth=1)
#inset_residuals.set_title('Residuals', fontsize=9)
inset_residuals.set_xlabel('q**2', fontsize=7)
inset_residuals.set_ylabel('Residuals/sigma', fontsize=7)
inset_residuals.ticklabel_format(axis="x", style="sci", scilimits=(0,0), useMathText=True) #sci test
inset_residuals.xaxis.get_offset_text().set_fontsize(7)
inset_residuals.tick_params(axis='both', which='major', labelsize=7)
#inset_residuals.grid(True)


# C- lot Kratky curve (top-right)
axs[0, 1].plot(xkrat, ykrat, color=col3)
axs[0, 1].set_xlabel('(qRg)^2')
axs[0, 1].set_ylabel('(qRg)^2 * I(q) / I(0)')
axs[0, 1].set_title('Normalized Kratky Plot')
axs[0, 1].axvline(x=1.73, color='grey', linestyle='--')
axs[0, 1].axhline(y=1.1, color='grey', linestyle='--')
axs[0, 1].grid(True)

# D- Plot the Volume of Correlation (VC) and Cumulative Integral
ax_vc = axs[1, 1]  # Left y-axis for Volume of Correlation
ax_integral = ax_vc.twinx()  # Create a twin axes for the right y-axis

# Plot the Volume of Correlation
ax_vc.plot(q_full, yvc_plot, label='Volume of Correlation', color=col1)
ax_vc.set_xlabel('q')
ax_vc.set_ylabel('(qRg)^2 * I(q) / I(0)', color=col1)
ax_vc.tick_params(axis='y', labelcolor=col1)
ax_vc.ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)#sci test
ax_vc.set_title('Volume of Correlation and Cumulative Integral')
ax_vc.axvline(x=0.3,color=col1, linestyle='--',alpha=0.8)
ax_vc.axvline(x=8/Rg,color=col2, linestyle='--',alpha=0.8)
# Add annotation for molecular weight (MW) in the VC plot

# Define text parts with different colors and formatted with a space for thousands
vc_text1 = f'MW Vc(q<0.3) = {MW1:,.0f}'.replace(',', ' ')  # Blue, formatted with space
vc_text2 = f'MW Vc(q<8/Rg) = {MW_alt:,.0f}'.replace(',', ' ')  # Red, formatted with space
vc_text3 = f'MW Porod (q<8/Rg) = {MV_P:,.0f}'.replace(',', ' ')  # Red, formatted with space
# Create a background box for the text
bbox_props = dict(boxstyle='round,pad=0.3', facecolor=col5, alpha=0.5)

# Add text with colored box for Volume of Correlation
ax_vc.text(0.17, 0.29, vc_text1,
           transform=ax_vc.transAxes, fontsize=10,
           verticalalignment='top', color=col1,
           bbox=bbox_props)  # Color for first mass

ax_vc.text(0.17, 0.22, vc_text2,
           transform=ax_vc.transAxes, fontsize=10,
           verticalalignment='top', color=col2,
           bbox=bbox_props)  # Color for second mass
ax_vc.text(0.17, 0.15, vc_text3,
           transform=ax_vc.transAxes, fontsize=10,
           verticalalignment='top', color=col4,
           bbox=bbox_props)  # Color for Porod

# Plot the cumulative integral on the right y-axis
ax_integral.plot(q_full, q_integral_cumulative, color=col3, label='Cumulative Integral')
ax_integral.ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)#sci test
ax_integral.set_ylabel('Integral of qI(q) dq', color=col3)
ax_integral.tick_params(axis='y', labelcolor=col3)

# Adding legends for both axes
ax_vc.legend(loc='lower left')
ax_integral.legend(loc='lower right')
ax_vc.grid(True)

###############################################
# Adjust layout to avoid overlap
#plt.tight_layout()
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05) #marge on the figure  
plt.savefig(output_file6)  # Save the graph as a PNG image
plt.savefig(output_file7,format="svg")  # Save the graph as a svg image
# Show the figure
plt.show()