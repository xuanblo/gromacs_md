rm -f _gmxtools*

echo -e "4 4" | gmx rms -s md_0_1.tpr -f md_0_1_whole_nojump.xtc -o rmsd.xvg -tu ns

# Calculate RMSD when not in equilibrium
echo -e "4 4" | gmx rms -s em.tpr -f md_0_1_whole_nojump.xtc -o rmsd_xtal.xvg -tu ns

# Analyze the radius of gyration of the simulated lysozyme using gyrate:
echo -e "1" | gmx gyrate -s md_0_1.tpr -f md_0_1_whole_nojump.xtc -o gyrate.xvg

# Grouping is required because only hydrogen bonds between chain a and chain b are counted
# The following is the method for creating an index

# https://jerkwin.github.io/2021/06/19/GROMACS%E6%B0%A2%E9%94%AE%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7hbond%E7%9A%84%E4%BD%BF%E7%94%A8%E5%8F%8A%E6%89%A9%E5%B1%95/
# 6m0j
#gmx make_ndx -f md_0_1.gro -o index.ndx
#> ri 1-597
#> name 18 Ace2
#> ri 598-791
#> name 19 S

#echo -e "ri 1-227\nname 18 ace2\nri 228-421\nname 19 S\nq" | gmx make_ndx -f md_0_1.gro -o index.ndx
receptor_pos=$1
receptor_name=$2
ligand_pos=$3
ligand_name=$4
cpu=$5

echo -e "ri $receptor_pos\nname 18 $receptor_name\nri $ligand_pos\nname 19 $ligand_name\nq" | gmx make_ndx -f md_0_1.gro -o index.ndx

# Obtain and analyze the root mean square fluctuation (RMSF) of the alpha-C atoms of the entire system
# RMSF calculates the fluctuation of each atom relative to its average position, representing the average change of the structure over time, and provides a characterization of the flexibility of different regions of the protein, corresponding to the temperature factor (B-factor) in crystallography. Typically, we expect RMSF to be similar to the temperature factor, which can be used to examine whether the simulation results are consistent with the crystal structure.
echo -e "19" | gmx rmsf -s md_0_1.tpr -f md_0_1_whole_nojump.xtc -o rmsf.xvg -res -n index.ndx

# Hydrogen bond count protein vs protein
echo -e "18 19" | gmx hbond -s md_0_1.tpr -f md_0_1_whole_nojump.xtc -num hbond_num.xvg -hbn -hbm -n index.ndx -nthreads $cpu -tu ns

bash ../../script/gmx_hbdat.bsh -s md_0_1.tpr -n hbond.ndx -m hbmap.xpm

# Calculate free energy
bash ../../script/gmx_mmpbsa.bsh -f md_0_1_whole_nojump.xtc -s md_0_1.tpr -n index.ndx -com Protein -pro $receptor_name -lig $ligand_name
