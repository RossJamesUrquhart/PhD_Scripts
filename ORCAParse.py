# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 17:42:17 2022

@author: avtei
"""
import numpy as np
import os, glob, sys
import matplotlib.pyplot as plt
from ase.io import read


# 1 Bohr = 0.52917724900001 Angstrom

def readin(fname):
    f = open(fname)
    content = f.read()
    f.close()
    return content


def fit_rms(ref_c,c):
    # move geometric center to the origin
    ref_trans = np.average(ref_c, axis=0)
    ref_c = ref_c - ref_trans
    c_trans = np.average(c, axis=0)
    c = c - c_trans

    # covariance matrix
    C = np.dot(c.T, ref_c)

    # Singular Value Decomposition
    (r1, s, r2) = np.linalg.svd(C)

    # compute sign (remove mirroring)
    if np.linalg.det(C) < 0:
        r2[2,:] *= -1.0
    U = np.dot(r1, r2)
    return (c_trans, U, ref_trans)

def calc_rmsd(c1, c2):
    rmsd = 0.0
    c_trans, U, ref_trans = fit_rms(c1, c2)
    new_c2 = np.dot(c2 - c_trans, U) + ref_trans
    rmsd = np.sqrt( np.average( np.sum( ( c1 - new_c2 )**2, axis=1 ) ) )
    return rmsd


class ORCAParse:
    def ValidateOutput(self):
        if not "***ORCA TERMINATED NORMALLY***" in self.raw:
            print(self.fname, "orca did not terminate normally!")
            self.valid = False
        elif "The optimization did not converge but reached the maximum number of" in self.raw:
            print(self.fname, "hit geom MaxIter!")
            self.valid = False
        else:
            self.valid = True

    def thermodynamics(self):
        #Search for the INNER ENERGY section
        pass
    
    def parse_energies(self):
        self.energies = np.ndarray((0,), np.float64)
        for part in self.raw.split("FINAL SINGLE POINT ENERGY")[1:]:
            part = float(part.split("\n")[0].strip())
            self.energies = np.hstack((self.energies, [part]))
        self.r_energies = self.energies - self.energies.min()
    def parse_dispersion(self):
        splits = self.raw.split("Dispersion correction")[1:]
        self.dispersions = np.ndarray((len(splits),))
        for i in range(len(splits)):
            E_disp = splits[i].split("\n")[0].strip()
            self.dispersions[i] = float(E_disp)
        
    def get_coords(self):
        self.coords = []
        self.atoms = []
        frames = self.raw.split("CARTESIAN COORDINATES (ANGSTROEM)")[1:]
        for i,frame in enumerate(frames):
            positions = []
            frame = frame.split("CARTESIAN COORDINATES (A.U.)")[0]
            for line in frame.split("\n"):
                line = line.split()
                if len(line) != 4:
                    continue
                positions.append([float(x) for x in line[1:]])
                if i == 0:
                    self.atoms.append(line[0])
            if i == 0:
                self.coords = np.array(positions).reshape(1, -1, 3)
            else:
                self.coords = np.vstack((self.coords, np.array(positions).reshape(1,-1,3)))
    
    def scan_bond(self, a0, a1):
        distances = np.ndarray((self.coords.shape[0]))
        for i in range(self.coords.shape[0]):
            positions = self.coords[i]
            d = np.linalg.norm(positions[a0] - positions[a1])
            distances[i] = d
        return distances
            
    def get_freqs(self):
        self.frequencies = []
        frames = self.raw.split("VIBRATIONAL FREQUENCIES")[1:]
        for i,frame in enumerate(frames):
            frequencies = []
            frame = frame.split("NORMAL MODES")[0]
            for line in frame.split("\n"):
                if "cm" not in line:
                    continue
                line = line.split()
                freq = float(line[1])
                frequencies.append(freq)
            if i == 0:
                self.frequencies = np.array(frequencies).reshape(1, -1)
            else:
                self.frequencies = np.vstack((self.frequencies, np.array(frequencies).reshape(1,-1)))
    def parse_free_energy(self):
        self.Gibbs = float(self.raw.split("Final Gibbs free energy         ...")[1].split("\n")[0].strip().split()[0])
        
    def time(self):
        time_str = self.raw.split("TOTAL RUN TIME:")[1].strip().split()
        days = int(time_str[0])
        hours = int(time_str[2])
        minutes = int(time_str[4])
        seconds = int(time_str[6])
        miliseconds = int(time_str[8])
        
        hours = hours + (days*24)
        minutes = minutes + (hours*60)
        seconds = seconds + (minutes*60)
        seconds = seconds + (miliseconds/1000)
        return seconds
    
    def __init__(self, fname):
        self.fname = fname
        self.raw = readin(fname)
        self.ValidateOutput()
        


    
if __name__ == "__main__":
    #fname = "K-DHP/TS/Insertion/ORCA_SCAN/Insertion_TS.out"
    #fname = "K-DHP/TS/Insertion/ORCA_TS_fromDNN/DNN_Berny.out"
    #fname = "K-DHP/TS/Deprotonation/ORCA_SCAN/Deprotonation_ScanTS.out"
    
    files = [x for x in glob.glob("Validation/Scan/IrClMe2_angle/*.out") if "atom77" not in x]
    for TS_out in files:
        print(TS_out)

        op = ORCAParse(TS_out)
        op.get_coords()
        print(op.coords.shape)
        op.parse_energies()
        print("op.valid:", op.valid)
        #print("op.valid:", op.energies)
        #plt.plot(op.r_energies, label=os.path.basename(TS_out))
        #plt.show()

        
        #op.get_freqs()
        #freqs = op.frequencies[-1]
        #print("Negative Freqs:", np.where(freqs < 0)[0].shape[0])
        
        a="""
        dists = op.scan_bond(a0=15, a1=9)
        
        plt.plot(dists, op.r_energies)
        plt.xlabel("Bond scan distance")
        plt.ylabel("rEnergy")
        #"""
        

        
