#!bin/bash/python

# Python script for running a benchmark for the Huang Coarse-Grained Model of Poly-3-Hexyl Thiophene on comet


import numpy as np
import os
import sys


import PolyModelFunctions_C60 as poly
import Huangparameters as Huang
import CG_System
#import RunLammps






def main():
    
    Script, deg_polymerization, num_polymer, box_length, zbox_length, dispersity, substrate, Solv_Screen, Ramp_Time, Hold_Time = sys.argv
    deg_polymerization = int(deg_polymerization)
    num_polymer = int(num_polymer)
    box_length = float(box_length)
    dispersity = float(dispersity)
    substrate = bool(substrate)
    Solv_Screen = float(Solv_Screen)
    zbox_length = float(zbox_length)
	
    Position_M, Position_S1, Position_S2, Bases, zBoxLength  =  poly.Gen_Many_Polymers( deg_polymerization, num_polymer, Huang.BondLengths, Huang.Angles, Huang.SigmaM_M, box_length, zbox_length, dispersity, substrate)
    
    Name = 'P3HT_Solvent_%d_%d' % ( deg_polymerization, num_polymer)
    Data_Filename = 'data.P3HT_Solvent_%d_%d' % ( deg_polymerization, num_polymer)


    CG_System.P3HTHuangDataLammps(Position_M, Position_S1, Position_S2, Data_Filename, box_length, zBoxLength)
    CG_System.Run_Huang_Equil(Name, Data_Filename, Solv_Screen, Ramp_Time = Ramp_Time, Hold_Time = Hold_Time)

    return

if __name__=='__main__': main() 
