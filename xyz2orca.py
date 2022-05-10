# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:30:09 2021

@author: avtei
"""
import os, sys, datetime, glob, time
from ase.io import read, write
import numpy as np

def buffer(string, L, end=True):
    while len(string) < L:
        string = string+"0"
    return string


CPU = """#!/bin/bash
#SBATCH --export=ALL
#SBATCH --output={}.log
#SBATCH --job-name={}
#SBATCH --account=tuttle-rmss
#SBATCH --partition=standard
#SBATCH --time="{}"
#SBATCH --ntasks={} --nodes=1
module purge
module load orca/5.0.1
module load orca/5.0.2

/opt/software/scripts/job_prologue.sh 

"""

def xyz2orca(fname):
    frames = read(fname, index=":")
    if len(frames) == 0:
        return False
    name = fname.split("/")[-1].replace(".xyz", "")
    species = frames[0].get_chemical_symbols()
    print("".join(species))
    folder = os.path.dirname(fname)
    
    CPUS = 4
    
    strT = '240:00:00'
    print(strT)
    
    for i,frame in enumerate(frames):
        jobname = fname.replace(".xyz", ".inp")
        inp = jobname.split("/")[-1]
        out = jobname.split("/")[-1].replace(".inp", ".out")
        if os.path.exists(folder+"/"+inp):
            print(f"Found: {out}, skipping.")
            return 0
        sbatch = open(jobname.replace(".inp", ".sh"), 'w')
        sbatch.write(CPU.format(jobname.replace(".inp", ""), jobname.replace(".inp", ""), strT, str(CPUS)))
        
        f = open(jobname, 'w')
        print(jobname)
        f.write(f"""# ORCAForces - {jobname}
# Basic Mode
#
! OPT wB97X D4 Def2-tzvpp  SlowConv TightSCF DEFGRID2 Def2/J RIJCOSX
%maxcore 3375
%PAL NPROCS {CPUS} END
%geom 
ReducePrint false
maxiter 1
END
* xyz 0 1
""".format(name))
        
        for coord, atom in zip(frame.get_positions(), species):
            line = [" "]*51
            line[0] = atom
            line[11:19] = buffer(str(coord[0]), 8)
            line[27:35] = buffer(str(coord[1]), 8)
            line[43:51] = buffer(str(coord[2]), 8)
            f.write("".join(line))
            f.write("\n")
        f.write("*") 
        f.close()
        
        sbatch.write("/opt/software/orca/5.0.2/orca {} > {}\n".format(inp, out))
        break
    
    sbatch.write("\n/opt/software/scripts/job_epilogue.sh\n")
    sbatch.close()
    
    cdir = os.path.abspath(".")
    os.chdir(folder)
    cmd = "sbatch {}".format(os.path.basename(jobname).replace(".inp", ".sh"))
    print(cmd)
    os.system(cmd)
    os.chdir(cdir)
    return True

    
def FilterChemistry(atomlist):
    species_order = ["H", "C", "N", "O", "F",  "P", "S", "Cl", "Ir"] # MUST BE ORDERED BY ATOMIC NUMBER
    for atom in atomlist:
        if atom not in species_order:
            return False
    return True
    
    

if __name__ == "__main__":
    tasks = glob.glob("PTM_structures/*.xyz")
    ordered_tasks = []
    ordered_tasks_natoms = []
    tasks_symbols = {}
    
    for task in tasks:
        try:
            mol = read(task)
        except:
            continue
        ordered_tasks.append(task)
        ordered_tasks_natoms.append(len(mol.get_chemical_symbols()))
        tasks_symbols[task] = mol.get_chemical_symbols()
    
    ordered_tasks = np.array(ordered_tasks)
    ordered_tasks_natoms = np.array(ordered_tasks_natoms)
    order_small_to_large = np.argsort(ordered_tasks_natoms)
    for task in ordered_tasks[order_small_to_large]:
        if FilterChemistry(tasks_symbols[task]):
            print(task, len(tasks_symbols[task]))
            xyz2orca(task)
        else:
            print(task, len(tasks_symbols[task]), "skipping because it contains a bad atom")
