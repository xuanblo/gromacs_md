rm -f *\# *.top *.tpr *.gro *.log *.cpt *.xtc *.edr *.itp *.xvg
rm -f *.trr
rm -f \_*

# General steps for any simulation include: building coordinates and topology, performing energy minimization, equilibration, running production MD simulation, and analyzing data collected during the simulation. We will follow these steps in this tutorial as well.

# 1 Prepare structure file, pdb file
# Remove small molecules to avoid errors
grep -v CL $1 | grep -v NAG | grep -v ZN > clean.pdb

# 2 Prepare topology file
# The processed PDB file contains only protein atoms without crystal water and other small molecules. We have also confirmed that no essential atoms are missing. Now we can generate the topology structure for it.
# top: topology, tip: include topology
# -ignh: ignore hydrogen atoms
# -ff: choose force field
gmx pdb2gmx -f clean.pdb -o processed.gro -ignh -ff amber99sb-ildn -water tip3p -missing
# OPLS-AA/L

# 3 Define the box and solvate
# Simulating a protein in vacuum is usually not biologically relevant. Simulating in a condensed phase is closer to reality. Therefore, the next step is to define a certain space around the protein and fill it with water.
# Define the box
# Add three different types of boxes and compare the number of water molecules in each box
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic
# Solvate with water
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top

# 4 Add ions
# Now grompp will integrate the coordinates, topology information, and parameters set in the .mdp file to generate a .tpr file.
cp ../../script/*.mdp ./
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
# We need to add counter ions to the system to ensure the total charge of the system is zero.
# Select solvent molecules for replacement
#gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral # 13
echo -e "SOL" | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral

# 5 Energy minimization
# Relax the solvated system to a lower energy state. This process is called energy minimization.
# Create input file for energy minimization
# The simulation box is too small
gmx grompp -debug -v -f minim.mdp -c ions.gro -p topol.top -o em.tpr -maxwarn 1
# Energy minimization (requires computational resources and takes a long time)
# 10^5-10^6 negative numbers
gmx mdrun -v -deffnm em -nt $2

# 5.1 Extract the relationship between the potential energy and the energy minimization steps
# View energy
echo -e "10 0" | gmx energy -f em.edr -o potential.xvg

# 6 NVT equilibration
# To start the actual dynamics simulation, we need to equilibrate the solvent and ions around the protein.
# Temperature and pressure
# Particle number (N), volume (V), average temperature (T), average pressure (P)
# Perform 100 ps NVT equilibration, approximately less than 1 hour with 16 cores
# First, use grompp to create the .tpr file for NVT equilibration
# Warm up from 0K to 300K
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
# -nt: number of threads, -nb: cpu/gpu
gmx mdrun -deffnm nvt -nt $2
# View energy
echo -e "16 0" | gmx energy -f nvt.edr -o hot_temperature.xvg

# 7 NPT equilibration
# First, use grompp to create the .tpr file for NPT equilibration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -nt $2

# 7.1 Analyze the pressure changes using the energy module to extract the pressure time series
echo -e "18 0" | gmx energy -f npt.edr -o pressure.xvg
# 7.2 Use the energy module again to extract the density time series from the npt.edr file
echo -e "24 0" | gmx energy -f npt.edr -o density.xvg

# 8 Production MD simulation
# Run grompp to generate the .tpr input file. We need to provide the coordinates and checkpoint file (which contains pressure coupling information) from NPT equilibration to ensure an accurate continuation of the dynamics process.
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
# Run
gmx mdrun -deffnm md_0_1 -nt $2

# Maintain integrity
echo -e "0" | gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_whole.xtc -pbc whole
echo -e "0" | gmx trjconv -s md_0_1.tpr -f md_0_1_whole.xtc -o md_0_1_whole_nojump.xtc -pbc nojump #-dt 50
