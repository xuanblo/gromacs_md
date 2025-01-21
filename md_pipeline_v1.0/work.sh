# -----------------
# step0 prepare   #
# -----------------
# create a folder with a given pdb file directory. per file per folder
python ../script/pdb2folder.py ./


# -----------------
# step1，run MD   #
# -----------------
nohup sh ../../script/step1.sh xxx.pdb 40


# --------------------------------------------
# step2，calculate RMSD, RMSF, HBOND, MMGBSA #
# --------------------------------------------
# check the order of the receptor and ligand
python ../../script/split_md_0_1_gro.py md_0_1.gro|less

# create index file and run gmx
sh ../../script/step2.sh 1-597 ace2 598-791 s 40
# note，如果有插入或者缺失的话，需要修改位置
#sh ../../script/step2.sh 1-597 ace2 598-790 s 40


# ------------------------------
# step3，statistic and plot    #
# ------------------------------
python ../script/merge_traj_rmsd_plot.py ./ rmsd.xlsx rmsd.pdf
python ../script/merge_res_rmsf_plot.py ./ rmsf.xlsx rmsf.pdf
python ../script/merge_hbond_plot.py ./ hbond.xlsx hbond.pdf
python ../script/mmgpbs_plot.py ./ traj_res_energe.xlsx mmpbsa_barplot.pdf


# -----------------------------------------------------------------
# step4，Calculate the binding energy of individual residues      #
# -----------------------------------------------------------------
python ../script/merge_traj_res_mmpbsa.py ./ merged_mmpbsa.csv merged_res.csv merge_traj_res_mmpbsa.pdf


# -----------------
# additional      #
# -----------------

# Extract trajectory
# Every 2fs per frame, 10ns, 30ns, 50ns
# 10ns
receptor="Sotrovimab_RBD"
for i in `ls -d */|cut -d '_' -f1`;do echo "echo" -e \"1\"\|gmx trjconv -s $i\_$receptor/md_0_1.tpr -f $i\_$receptor/md_0_1_whole_nojump.xtc -o $i\_$receptor/$i\_md_0_1.10.pdb -dump 10000;done >run_trjconv.sh
# 30ns
for i in `ls -d */|cut -d '_' -f1`;do echo "echo" -e \"1\"\|gmx trjconv -s $i\_$receptor/md_0_1.tpr -f $i\_$receptor/md_0_1_whole_nojump.xtc -o $i\_$receptor/$i\_md_0_1.30.pdb -dump 30000;done >>run_trjconv.sh
# 50ns
for i in `ls -d */|cut -d '_' -f1`;do echo "echo" -e \"1\"\|gmx trjconv -s $i\_$receptor/md_0_1.tpr -f $i\_$receptor/md_0_1_whole_nojump.xtc -o $i\_$receptor/$i\_md_0_1.50.pdb -dump 50000;done >>run_trjconv.sh



# Analyze salt bridges between residues in the simulation. The program will output a series of xvg files, providing distances between -/-, +/- (most relevant), and +/+ residues.
gmx saltbr -s md_0_1.tpr -f md_0_1_whole_nojump.xtc -dt 10000


# Calculate hydrogen bond occupancy
for i in `ls -d */|cut -d '_' -f1`;do python ../script/resn2hbond.py $i\_ACE2_S1/md_0_1.gro $i\_ACE2_S1/hbdat.dat $i;done

