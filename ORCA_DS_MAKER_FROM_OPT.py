#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 00:38:44 2021

@author: rkb19187
"""

import h5py, glob, os, sys
import re
import numpy as np
from ORCAParse import *
from torchani.utils import ChemicalSymbolsToInts

# 1 Bohr = 0.52917724900001 Angstrom

protons = {1: "H",
           5: "B",
           6: "C",
           7: "N",
           8: "O",
           9: "F",
           15:"P",
           16:"S",
           17:"Cl",
           77: "Ir"}

def readin(fname):
    f = open(fname, 'r')
    content = f.read()
    return content

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def atom_list_out(outfile):
    # Get list of atoms
    outcontent = readin(outfile)
    atomlist = []
    outcontent = outcontent.split("CARTESIAN COORDINATES (ANGSTROEM)")[1]
    outcontent = outcontent.split("CARTESIAN COORDINATES (A.U.)")[0]
    for line in outcontent.split("\n"):
        line = line.split()
        if len(line) != 4:
            continue
        atomlist.append(line[0])
    return atomlist

def get_opt(optfile):
    content = readin(optfile)
    #Validate the content
    if not "$coordinates" in content:
        print("Missing $coordinates")
        return [False]*5
    elif not "$energies" in content:
        print("Missing $energies")
        return [False]*5
    elif not "$gradients" in content:
        print("Missing $gradients")
        return [False]*5
    
    coord_data = content.split("$coordinates\n")[1].split("#")[0]
    line0 = coord_data.split("\n")[0]
    frames, atoms = line0.strip().split()
    frames, atoms = int(frames), int(int(atoms)/3)
    coord_data = " ".join(coord_data.split("\n")[1:]) #remove first line & and remove all newlines to make it one long line
    coord_data = coord_data.strip().split() # make into list of texts
    coord_data = np.array(coord_data) # convert list to numpy array
    coord_data = coord_data.astype(np.float64) # convert texts to numbers
    coordinates = np.ndarray((frames, atoms, 3))
    for i,frame in enumerate(np.split(coord_data, frames)): # Split it into the known number of frames we have
        # each frame has the shape atoms*3, we want to convert it to (atoms, 3)
        frame = np.vstack(np.split(frame, atoms))
        coordinates[i] = frame
    
    energy_data = content.split("$energies\n")[1].split("\n\n")[0]
    energy_data = " ".join(energy_data.split("\n")[1:]) #remove first line & and remove all newlines to make it one long line
    energy_data = energy_data.strip().split() # make into list of texts
    energy_data = np.array(energy_data) # convert list to numpy array
    energies = energy_data.astype(np.float64) # convert texts to numbers

    gradients_data = content.split("$gradients\n")[1].split("$")[0]
    gradients_data = " ".join(gradients_data.split("\n")[1:]) #remove first line & and remove all newlines to make it one long line
    gradients_data = gradients_data.strip().split() # make into list of texts
    gradients_data = np.array(gradients_data) # convert list to numpy array
    gradients_data = gradients_data.astype(np.float64) # convert texts to numbers
    forces = np.ndarray((frames, atoms, 3))    
    for i,frame in enumerate(np.split(gradients_data, frames)): # Split it into the known number of frames we have
        # each frame has the shape atoms*3, we want to convert it to (atoms, 3)
        frame = np.vstack(np.split(frame, atoms))
        forces[i] = frame

    return frames, atoms, energies, coordinates, forces


species_order = ["H", "B", "C", "N", "O", "F", "P", "S", "Cl", "Ir"] # MUST HAVE ONE OF EACH
check = {}
for s in species_order:
    check[s] = False
    
if __name__ == "__main__":
    species_to_indices = ChemicalSymbolsToInts(species_order)
    #user = "Alex" if os.environ['COMPUTERNAME'] ==  'PAC-FOG-PC-0028' else "Ross"
    user = "Alex"
    print(user)
    
    opts = glob.glob("DatasetFolder/Jan2022Data/*.opt")
    np.random.shuffle(opts)
    
    mol, E, C, S, F = [],[],[],[],[]
    
    hdf5file = "N_vs_Error/Unknown/N=Unknown.h5"
    outfolder = os.path.dirname(hdf5file)
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
    master = h5py.File(hdf5file, 'w')
    
    N = 0
    target = 1000
    for optfile in opts:
        outfile = optfile.replace(".opt", ".out")
        if not os.path.exists(outfile):
            print(f"Couldn't find .out file: {outfile}")
            continue
        
        atomlist = atom_list_out(outfile)
        badatom = False
        for a in np.unique(atomlist):
            if a not in species_order:
                badatom = True
                break
        if "I" in atomlist:
            sys.exit()
        if badatom:
            print("Atom that doesnt exist in species_order in:", atomlist)
            continue
            
        print(optfile, end=" ")
        frames, atoms, energies, coordinates, forces = get_opt(optfile)
        if frames == False:
            print("Bad file")
            continue
        
                
        print(atomlist)
        #index_tensor = species_to_indices(atomlist)
        #print(index_tensor)
        
        N += forces.shape[0]    
        print(coordinates.shape[0])

        name = os.path.basename(optfile.replace(".opt", ""))
        L = coordinates.shape[0]

        mol.append(master.create_group(name))
        
        species = np.array(atomlist, dtype = h5py.special_dtype(vlen=str) ) 
        print(species)
        coordinates = coordinates * 0.52917724900001
        forces = forces * 0.52917724900001
        print(name)
    
        E.append(mol[-1].create_dataset("energies", (energies.shape[0],), dtype='float64'))
        E[-1][()] = energies
        
        C.append(mol[-1].create_dataset("coordinates", coordinates.shape, dtype='float32'))
        C[-1][()] = coordinates
        
        F.append(mol[-1].create_dataset("forces", forces.shape, dtype='float32'))
        F[-1][()] = forces
        
        S.append(mol[-1].create_dataset("species", data=species))
        #S.append(mol[-1].create_dataset("species", data=index_tensor.numpy()))
        
        if N >= target:
            print(f"N = {N}, which is >= {target}, stopping...")
            break


    master.close()

if not os.path.exists(f"N_vs_Error/{N}"):
    os.mkdir(f"N_vs_Error/{N}")
os.rename(hdf5file, hdf5file.replace("Unknown", str(N)))
os.rmdir(f"N_vs_Error/Unknown")