import glob
import os
import sys
import shutil

pdb_dir = sys.argv[1]

pdbs = glob.glob(os.path.join(pdb_dir, '*.pdb'))

for pdb in pdbs:
    pdf_name = os.path.basename(pdb)
    dir_name = pdf_name.replace('.pdb', '')
    dir_path = os.path.join(pdb_dir, dir_name)
    os.makedirs(dir_path, exist_ok=True)
    shutil.move(pdb, dir_path)