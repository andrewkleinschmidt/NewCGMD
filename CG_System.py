# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 16:44:34 2015
File containing a function for preparing a LAMMPS Data file to implement the huang model
@author: Samuel
"""

import sys
import numpy as np
from Huangparameters import *
import os
import matplotlib.pyplot as plt
import subprocess
import Configure
import time
import math

def P3HTHuangDataLammps(Position_M, Position_S1, Position_S2, Filename, box_length, zBoxLength = Box_Length):

    Atom_Masses = {}
    Bond_Param = {}
    Angle_Param = {}
    Dihed_Param = {}
    Improp_Param = {}
    All_Keys = []
    
    with open(Configure.Param_Path + "P3HT.txt", 'r') as f:
        Divisions_Temp = f.readline().split(',')
        Divisions = []
        for d in Divisions_Temp:
            Divisions.append(int(d.strip()))
        NumAtomTypes = Divisions[0]
        NumBondTypes = Divisions[1]
        NumAngleTypes = Divisions[2]
        NumDihedralTypes = Divisions[3]
        NumImproperTypes = Divisions[4]
        for i in range(NumAtomTypes):
            key = f.readline().strip()
            val = float(f.readline().strip())
            Atom_Masses[key] = val
            All_Keys.append(key)
        for i in range(NumBondTypes):
            params = []
            key = f.readline().strip()
            params_temp = f.readline().split(',')
            for p in params_temp:
                params.append(float(p.strip()))
            Bond_Param[key] = params
                             #All_Keys.append(key)
        for i in range(NumAngleTypes):
            params = []
            key = f.readline().strip()
            params_temp = f.readline().split(',')
            for p in params_temp:
                params.append(float(p.strip()))
            Angle_Param[key] = params
                             #All_Keys.append(key)
        for i in range(NumDihedralTypes):
            params = []
            key = f.readline().strip()
            params_temp = f.readline().split(',')
            for p in params_temp:
                params.append(float(p.strip()))
            Dihed_Param[key] = params
            All_Keys.append(key)
        for i in range(NumImproperTypes):
            params = []
            key = f.readline().strip()
            params_temp = f.readline().split(',')
            for p in params_temp:
                params.append(float(p.strip()))
            Improp_Param[key] = params
                             #All_Keys.append(key)
                
    File = open(Filename, 'w')
    
    NumChains = np.size(Position_M, 0)
    #NumAtoms = NumChains*ChainLength*3 + NumPCBM
    #NumBonds = NumChains*(ChainLength - 1 + ChainLength*2)
    #NumAngles = NumChains*(ChainLength + 3*ChainLength - 4)
    #NumDih = NumChains*(ChainLength - 3 + ChainLength - 1 + ChainLength - 1 + ChainLength - 1)
    #NumImproper = NumChains*(ChainLength - 2)
    #NumImproper = 0
    NumAtoms = 0
    NumBonds = 0
    NumAngles = 0
    NumDih = 0
    NumImproper = 0
    for polymer in Position_M:
        NumAtoms += 3 * np.size(polymer,0)
        NumBonds += 3 * np.size(polymer,0) - 1
        NumAngles += 4 * np.size(polymer,0) - 4
        NumDih += 4 * np.size(polymer,0) - 6
        NumImproper += np.size(polymer,0) - 2

    File.write('Created by Sam Root\n\n')
    File.write('\t%d\tatoms\n' % NumAtoms)
    File.write('\t%d\tbonds\n' % NumBonds)
    File.write('\t%d\tangles\n' % NumAngles)
    File.write('\t%d\tdihedrals\n'% NumDih)
    File.write('\t%d\timpropers\n\n' % NumImproper)

    File.write('\t%d\tatom types\n' % NumAtomTypes)
    File.write('\t%d\tbond types\n' % NumBondTypes)
    File.write('\t%d\tangle types\n' % NumAngleTypes)
    File.write('\t%d\tdihedral types\n' % NumDihedralTypes)
    File.write('\t%d\timproper types\n\n' % NumImproperTypes)

    File.write('\t%.6f\t%.6f xlo xhi\n' %(0., box_length))
    File.write('\t%.6f\t%.6f ylo yhi\n' %(0., box_length))
    File.write('\t%.6f\t%.6f zlo zhi\n\n' %(0., zBoxLength))

    # Declare Masses
    File.write('Masses\n\n')
    for i in range(NumAtomTypes):
        File.write('\t%d %.4f\n' % (i+1, Atom_Masses[All_Keys[i]]))
    
    # Declare Dihedral Coeffs
    File.write('Dihedral Coeffs\n\n')
    for i in range(NumDihedralTypes):
        File.write('\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' % (i+1, Dihed_Param[All_Keys[i+NumAtomTypes]][0], Dihed_Param[All_Keys[i+NumAtomTypes]][1], Dihed_Param[All_Keys[i+NumAtomTypes]][2], Dihed_Param[All_Keys[i+NumAtomTypes]][3], Dihed_Param[All_Keys[i+NumAtomTypes]][4]))

    # Declare atoms initial conditions
    File.write('Atoms\n\n')
    Atom_id = 0
    #Mol_id = 0
    #k = 0
    """for i in range(NumChains):
        Mol_id += 1
        for j in range(ChainLength):
        Atom_id += 1
        File.write('%d %d 1 %.2f %.2f %.2f\n' % ( Atom_id, Mol_id, Position_M[k,0], Position_M[k,1], Position_M[k,2]))
        Atom_id += 1
        File.write('%d %d 2 %.2f %.2f %.2f\n' % ( Atom_id, Mol_id,  Position_S1[k,0], Position_S1[k,1], Position_S1[k,2]))
        Atom_id += 1
        File.write('%d %d 3 %.2f %.2f %.2f\n' % ( Atom_id, Mol_id,  Position_S2[k,0], Position_S2[k,1], Position_S2[k,2]))
        k += 1"""
    for Mol_id, polymer in enumerate(Position_M):
        for atom in range(np.size(polymer,0)):
            Atom_id += 1
            File.write("%d %d 1 %.2f %.2f %.2f\n" % (Atom_id, Mol_id+1, Position_M[Mol_id][atom][0], Position_M[Mol_id][atom][1], Position_M[Mol_id][atom][2]))
            Atom_id += 1
            File.write("%d %d 2 %.2f %.2f %.2f\n" % (Atom_id, Mol_id+1, Position_S1[Mol_id][atom][0], Position_S1[Mol_id][atom][1], Position_S1[Mol_id][atom][2]))
            Atom_id += 1
            File.write("%d %d 3 %.2f %.2f %.2f\n" % (Atom_id, Mol_id+1, Position_S2[Mol_id][atom][0], Position_S2[Mol_id][atom][1], Position_S2[Mol_id][atom][2]))
    """for i in range(NumPCBM):
        Mol_id += 1
        Atom_id +=1
        File.write('%d %d 4 %.2f %.2f %.2f\n' % (Atom_id, Mol_id, Position_C60[i,0], Position_C60[i,1], Position_C60[i,2]))"""


    # Declare Bonding Topology
    File.write('\n\nBonds\n\n')
    E = 0
    B =1
    for polymer in Position_M:
        for j in range(np.size(polymer,0) - 1):
            E+= 1
            File.write('%d 1 %d %d\n' % (E, B, B+3))
            E += 1
            File.write('%d 2 %d %d\n' % (E, B, B+1))
            E += 1
            File.write('%d 3 %d %d\n' % (E, B+1, B+2))
            B += 3
        E += 1
        File.write('%d 2 %d %d\n' % (E, B, B+1) )
        E += 1
        File.write('%d 3 %d %d\n' % (E, B+1, B+2))
        B+= 3


    # Declare Angle Topology
    File.write('\n\nAngles\n\n')
    E = 0
    B = 1
    for polymer in Position_M:
        for j in range(np.size(polymer,0) - 2):
            E += 1
            File.write('%d 1 %d %d %d\n' % (E, B, B+3, B+6))
            E +=1
            File.write('%d 2 %d %d %d\n' % (E, B, B+1, B+2))
            E +=1
            File.write('%d 3 %d %d %d\n' % (E, B, B+3, B+4))
            E += 1
            File.write('%d 4 %d %d %d\n' % (E, B+1, B, B+3))
            B+=3
        E += 1
        File.write('%d 2 %d %d %d\n' % (E, B, B+1, B+2))
        E += 1
        File.write('%d 3 %d %d %d\n' % (E, B, B+3, B+4))
        E += 1
        File.write('%d 4 %d %d %d\n' % (E, B+1, B, B+3))
        E += 1
        File.write('%d 2 %d %d %d\n' % (E, B+3, B+4, B+5))
        B+= 6

    # Declare Dihedral Topology
    File.write('\n\nDihedrals\n\n')
    E = 0
    B = 1
    for polymer in Position_M:
        for j in range(np.size(polymer,0) - 3):
            E += 1
            File.write('%d 1 %d %d %d %d\n' % (E, B, B+3, B+6, B+9) )
            E += 1
            File.write('%d 2 %d %d %d %d\n' % (E, B+1, B, B+3, B+4))
            E += 1
            File.write('%d 3 %d %d %d %d\n' % (E, B, B+3, B+4, B+5))
            E += 1
            File.write('%d 4 %d %d %d %d\n' % (E, B+2, B+1, B, B+3))
            B +=3
        E += 1
        File.write('%d 2 %d %d %d %d\n' % (E, B+1, B, B+3, B+4))
        E += 1
        File.write('%d 3 %d %d %d %d\n' % (E, B, B+3, B+4, B+5))
        E += 1
        File.write('%d 4 %d %d %d %d\n' % (E, B+2, B+1, B, B+3))
        B+= 3
        E+= 1
        File.write('%d 2 %d %d %d %d\n' % (E, B+1, B, B+3, B+4))
        E+= 1
        File.write('%d 3 %d %d %d %d\n' % (E, B, B+3, B+4, B+5))
        E += 1
        File.write('%d 4 %d %d %d %d\n' % (E, B+2, B+1, B, B+3))
        B += 6


    # Declare Improper Dihedral Topology
    E = 0
    B = 1
    File.write('\n\nImpropers\n\n')
    for polymer in Position_M:
        for j in range(np.size(polymer,0) -2):
            E += 1
            File.write('%d 1 %d %d %d %d\n' % (E, B, B+3, B+6, B+7))
            B += 3
        B += 6

    File.close()
        
    return


def Run_Huang_Equil(Name, Data_In, Solv_Screen, Dump = 'None', Num_Steps= 1000000, Time_Step = 4, Ramp_Time = 100, Hold_Time = 0, Solv_Screen_End = 1.0, NVT = True, Substrate = True):
    """
    Function for generating an input file for an NVT equilibration simulation of Huang Model
    
    Input Arguments:
        Required:
        In_File = name of the input file to be created
        Restart_Out = name of restart file to be outputed
        Data_File | Restart_In = File for declaring initial position and topology (atleast on of these is required)
        Pair_table = name of pair table file
    Optional Arguments:
        Dump = name of dumpfile to be outputed, Default is none
        Num_Steps = number of timesteps to integrate equations of motion
        Time_Step = size of time increment (fs)
        Ramp_Time = total time running while ramping solvent screening (ns)
        Hold_Time = total time running after ramping solvent screening (ns)
        Temp = Temperature to run simulation at
    """
    Init_Temp = Configure.Template_Path + "init_Langevin.in"
    Sub_Temp = Configure.Template_Path + "sub_Lammps"
    Sub_Temp = Configure.Template_Path + "GPU_sub"
    Init_File = "in.init_" + Name
    Sub_File = "sub_" + Name
    File_Out1 = "log.init_%s" % Name
    Restart_Out = "restart.%s_800_1" % Name
    Nodes = 4
    NProcs = 24
    if Substrate:
        Sub_Bound = 's'
    else:
        Sub_Bound = 'p'
    print "ssh",Configure.Comet_Login,"mkdir " + Configure.Comet_Path % Name

    cmd = "mkdir " + Configure.Comet_Path % Name
    subprocess.call(["ssh", Configure.Comet_Login, cmd])
    
    #Turn below into a function
    with open(Init_Temp, 'r') as f:
        template = f.read()
    s = template.format(Solv_Screen = Solv_Screen, NVT = NVT, Sub_Bound = Sub_Bound, Substrate = Substrate,Restart_Out = Restart_Out,Data_In = Data_In,Temp_In = 300, Temp_Out = 800)
    with open(Init_File, 'w') as f:
        f.write(s)
    """  code to determine # of processors, num nodes"""
    with open(Sub_Temp, 'r') as f:
        sub_template = f.read()
    s = sub_template.format(Sim_Name = "init_" + Name, path = Configure.Comet_Path % Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs) #EDIT
    with open(Sub_File, 'w') as f:
        f.write(s)
    os.system( Configure.c2c % (Sub_File, Name))
    os.system( Configure.c2c % (Init_File, Name))
    os.system( Configure.c2c % (Data_In, Name))
    os.system( Configure.c2l % (Name, File_Out1))
    try:
        File = open(File_Out1,'r')
    except:
        subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (Name, Sub_File)])

    Finished = False
    i = 0
    while not Finished:
        os.system( Configure.c2l % (Name, Restart_Out))
        try:
            File = open(Restart_Out,'r')
            Finished = True
        except:
            print "Sleeping process", i, "minutes"
            time.sleep(600)
            i += 10
    os.system( 'rm %s' % Restart_Out)

    Run_Temp = Configure.Template_Path + "Langevin_Ramp.in"
    if NVT:
        NVT_str = "NVT"
    else:
        NVT_str = "NPT"
    Stop = Num_Steps * int(math.ceil(Ramp_Time % (Num_Steps * Time_Step / 1000000)))

    for i in range (int(math.ceil(Ramp_Time % (Num_Steps * Time_Step / 1000000)))):
        Restart_In = Restart_Out
        Restart_Out = "restart.%s_%s_%d_Ramp" % Name, NVT_str, i
        Run_File = "in.%s_%s_%d_Ramp" % Name, NVT_str, i
        File_Out1 = "log.%s_%s_%d_Ramp" % Name, NVT_str, i
        with open(Run_Temp, 'r') as f:
            template = f.read()
        template.format(Solv_Screen_Start = Solv_Screen, Solv_Screen_End = Solv_Screen_End, NVT = NVT, Substrate = Substrate,Restart_Out = Restart_Out,Restart_In = Restart_In,Temp_In = 300, Temp_Out = 300, Time_Step = Time_Step, Num_Steps = Num_Steps, Stop = Stop)
        with open(Run_File, 'w') as f:
            f.write(template)
            """  code to determine # of processors, num nodes"""
        with open(Sub_Temp, 'r') as f:
            sub_template = f.read()
        sub_template.format(Sim_Name = NVT_str + "_" +  Name + "_d" % i, path = Configure.Comet_Path % Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
        with open(Sub_File, 'w') as f:
            f.write(sub_template)
        os.system( Configure.c2c % (Sub_File, Name))
        os.system( Configure.c2c % (Run_File, Name))
        os.system( Configure.c2c % (Restart_In, Name))
        os.system( Configure.c2l % (Name, File_Out1))
        try:
            File = open(File_Out1,'r')
        except:
            subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (Name, Sub_File)])
        Finished = False
        i = 0
        while not Finished:
            os.system( Configure.c2l % (Name, Restart_Out))
            try:
                File = open(Restart_Out,'r')
                Finished = True
            except:
                print "Sleeping process", i, "minutes"
                time.sleep(600)
                i += 10
        os.system( 'rm %s' % Restart_Out)

    Run_Temp = Configure.Template_Path + "Langevin_Hold.in"

    for i in range (int(math.ceil(Hold_Time % (Num_Steps * Time_Step / 1000000)))):
        Restart_In = Restart_Out
        Restart_Out = "restart.%s_%s_%d_Hold" % Name, NVT_str, i
        Run_File = "in.%s_%s_%d_Hold" % Name, NVT_str, i
        File_Out1 = "log.%s_%s_%d_Hold" % Name, NVT_str, i
        with open(Run_Temp, 'r') as f:
            template = f.read()
        template.format(Solv_Screen_End = Solv_Screen_End, NVT = NVT, Substrate = Substrate,Restart_Out = Restart_Out,Restart_In = Restart_In,Temp_In = 300, Temp_Out = 300, Time_Step = Time_Step, Num_Steps = Num_Steps)
        with open(Run_File, 'w') as f:
            f.write(template)
            """  code to determine # of processors, num nodes"""
    with open(Sub_Temp, 'r') as f:
        sub_template = f.read()
    sub_template.format(Sim_Name = NVT_str + "_" +  Name + "_d" % i, path = Configure.Comet_Path % Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
    with open(Sub_File, 'w') as f:
        f.write(sub_template)
    os.system( Configure.c2c % (Sub_File, Name))
    os.system( Configure.c2c % (Run_File, Name))
    os.system( Configure.c2c % (Restart_In, Name))
    os.system( Configure.c2l % (Name, File_Out1))
    try:
        File = open(File_Out1,'r')
    except:
        subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (Name, Sub_File)])
        Finished = False
        i = 0
        while not Finished:
            os.system( Configure.c2l % (Name, Restart_Out))
            try:
                File = open(Restart_Out,'r')
                Finished = True
            except:
                print "Sleeping process", i, "minutes"
                time.sleep(600)
                i += 10
            sub_template.format(Sim_Name = NVT_str + "_" +  Name + "_d" % i, path = Configure.Comet_Path % Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
            with open(Sub_File, 'w') as f:
                f.write(sub_template)
            os.system( Configure.c2c % (Sub_File, Name))
            os.system( Configure.c2c % (Run_File, Name))
            os.system( Configure.c2c % (Restart_In, Name))
            os.system( Configure.c2l % (Name, File_Out1))
            try:
                File = open(File_Out1,'r')
        except:
            subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (Name, Sub_File)])
            Finished = False
            i = 0
            while not Finished:
                os.system( Configure.c2l % (Name, Restart_Out))
                try:
                    File = open(Restart_Out,'r')
                    Finished = True
                except:
                    print "Sleeping process", i, "minutes"
                    time.sleep(600)
                    i += 10
    os.system( 'rm %s' % Restart_Out)

    """File.write( '# Input file for running an NVT Equilibration in LAMMPS, Filename = %s\n\n' % In_File)
    File.write('units real\n')
    File.write('atom_style molecular\n')
    File.write('bond_style table linear 200\n')
    File.write('pair_style lj/cut 25.0 \n')
    File.write('angle_style table linear 200\n')
    File.write('dihedral_style multi/harmonic\n')
    File.write('special_bonds lj 0 0 0\n')
    File.write('improper_style none\n')
    File.write('kspace_style none\n')
    File.write('read_data %s\n' % Data_File)
    #File.write('pair_coeff 1 1 HuangPair2.table P1P1\n')
    #File.write('pair_coeff 1 2 HuangPair2.table P1P2\n')
    #File.write('pair_coeff 1 3 HuangPair2.table P1P3\n')
    #File.write('pair_coeff 2 2 HuangPair2.table P2P2\n')
    #File.write('pair_coeff 2 3 HuangPair2.table P2P3\n')
    #File.write('pair_coeff 3 3 HuangPair2.table P3P3\n')
    #File.write('pair_coeff 1 4 HuangPair3.table P1F1\n')
    #File.write('pair_coeff 2 4 HuangPair3.table P2F1\n')
    #File.write('pair_coeff 3 4 HuangPair3.table P3F1\n')
    #File.write('pair_coeff 4 4 HuangPair3.table F1F1\n')
    File.write('bond_coeff 1 HuangBond.table Bond1\n' )
    File.write('bond_coeff 2 HuangBond.table Bond2\n' )
    File.write('bond_coeff 3 HuangBond.table Bond3\n' )
    File.write('angle_coeff 1 HuangAngle.table Angle1\n')
    File.write('angle_coeff 2 HuangAngle.table Angle2\n')
    File.write('angle_coeff 3 HuangAngle.table Angle3\n')
    File.write('angle_coeff 4 HuangAngle.table Angle4\n')
    File.write('pair_coeff 1 1 .4368 3.385\n')
    File.write('pair_coeff 2 2 .2948 4.517\n')
    File.write('pair_coeff 3 3 .03182 4.287\n')
    File.write('pair_coeff 4 4 3.234 9.355\n')
    File.write('pair_modify shift yes mix geometric\n')
    if (Dump != 'None'):
        File.write('dump 1 all custom 1000 %s id type mol xs ys zs vx vy vz\n' % Dump)
    File.write('neighbor 10.0 bin\n')
    File.write('neigh_modify every 10 delay 0 one 1000\n')
    File.write('velocity all create %f 1223\n' % Temp)
    File.write('thermo_style custom step temp press etotal density\n')
    File.write('fix 1 all nve/limit 0.5\n')
    File.write('fix 2 all langevin %f %f 15 1234\n' % (Temp, Temp));#  File.write('fix 3 all deform 1 x scale .08 y scale .08 z scale .08 remap x \n')
    if zBoxLength != Box_Length:
        File.write('fix 3 all wall/lj126 zlo EDGE .369 3.35 13.3047')
    File.write('timestep %d\n' % Time_Step)
    File.write('thermo 100\n')
    File.write('run %d\n' % Num_Steps)
    File.write('unfix 1\n')
    File.write('unfix 2\n')
    File.write('write_restart %s\n' % Restart_Out)
    File.close()"""
    
    #os.system('mpirun -np 8 lammps <  %s' % In_File )
    
    return
    
    
def EquilibrateNPT( In_File, Restart_In, Restart_Out, Dump = 'none', Num_Steps = 100000, Time_Step = 1, Temp = 423, Press = 1, NP = 8, Num_Threads=16):
    """ Function for equilibrating an MD Simulation in an NPT Ensemble

    Required Inputs:
            In_File = Name of simulation input file
            
            Restart_In = Name of Restart file to read in
            
            Restart_Out = Name of Restart file to write out
            
            
    Optional Inputs:
            Dump = name of dump file to output, Default = none
            
            Num_Steps = number of steps to integrate over
            
            Time_Step = time increment
            
            Temp = Simulation Temperature
            
            Press = Simulation Pressure
    """
    File = open( In_File, 'w')
    
    File.write('# Input file for running an NPT Equilibration with LAMMPS, File name = %s\n\n' % In_File)
    File.write('read_restart %s\n' % Restart_In)
    File.write('pair_style table linear 187 \n');File.write('angle_style table linear 200\n')
    File.write('bond_style table linear 200\n');File.write('special_bonds lj 0 0 0\n')
    File.write('pair_coeff 1 1 HuangPair2.table P1P1\n')
    File.write('pair_coeff 1 2 HuangPair2.table P1P2\n')
    File.write('pair_coeff 1 3 HuangPair2.table P1P3\n')
    File.write('pair_coeff 2 2 HuangPair2.table P2P2\n')
    File.write('pair_coeff 2 3 HuangPair2.table P2P3\n')
    File.write('pair_coeff 3 3 HuangPair2.table P3P3\n')
    File.write('pair_coeff 1 4 HuangPair3.table P1F1\n')
    File.write('pair_coeff 2 4 HuangPair3.table P2F1\n')
    File.write('pair_coeff 3 4 HuangPair3.table P3F1\n')
    File.write('pair_coeff 4 4 HuangPair3.table F1F1\n')
    File.write('bond_coeff 1 HuangBond.table Bond1\n' )
    File.write('bond_coeff 2 HuangBond.table Bond2\n' )
    File.write('bond_coeff 3 HuangBond.table Bond3\n' )
    File.write('angle_coeff 1 HuangAngle.table Angle1\n')
    File.write('angle_coeff 2 HuangAngle.table Angle2\n')
    File.write('angle_coeff 3 HuangAngle.table Angle3\n')
    File.write('angle_coeff 4 HuangAngle.table Angle4\n')
    #File.write('log log.NPT_%d_%d\n' % (NP, Num_Threads))
    File.write('neighbor 10.0 bin\n')
    File.write('neigh_modify every 10 delay 0 one 2000\n')
    File.write('fix 1 all npt temp %d %d 100 iso 0.0 %d 1000 drag 2\n' % (Temp, Temp, Press))
    File.write('fix 2 all momentum 1 linear 1 1 1\n')
    File.write('thermo_style custom step temp press etotal density\n')
    
    if (Dump != 'none'):
        File.write('dump 1 all custom 200 %s id type mol xs ys zs vx vy vz\n' % Dump)
    
    File.write('thermo 200\n')
    File.write('timestep %d\n' % Time_Step)
    File.write('run %d\n' % Num_Steps)
    File.write('unfix 1\n')
    File.write('unfix 2\n')
    File.write('write_restart %s\n' % Restart_Out )
    File.close()
    #os.system('export OMP_NUM_THREADS=%s' % Num_Threads)
    #os.system('mpirun -np %d lammps < %s ' % (NP, In_File))
    
    return

def Extract_Trajectory(Filename, N, T, K):
    """
    Function to extract the trajectories from a sorted lammpstrj file
    inputs: Filename
            N = number of particles of each type
            T = number of timesteps
            k = # of types of particles
            
    outputs numpy arrays
    """
    
    File = open(Filename)
    Positions = np.zeros((K, N*T, 3))
    i= 0
    j= 0
    k = 0
    for line in File:
        line = line.split()
        if len(line) == 9:
            if line[1] == '1':
                Positions[0, i] = [float(line[3]), float(line[4]), float(line[5])]
                i += 1
            if line[1] == '2':
                Positions[1, j] = [float(line[3]), float(line[4]), float(line[5])]
                j+=1
            if line[1] == '3':
                Positions[2, k] = [float(line[3]), float(line[4]), float(line[5])]
                k+=1
    return Positions
    
def Order_Parameter(Positions, N, T, k, chainlength, avgT):
    Num_Chains = N/chainlength
    Num_angles = (chainlength-2)*Num_Chains
    Orientationx = np.zeros((Num_angles, T))
    for i in range(T):
        a=0
        b=0
        for j in range(Num_Chains):
            for k in range(chainlength - 2):
                rjk= (Positions[0, a*i] - Positions[0, a*i+2])
                rjk = rjk/np.linalg.norm(rjk)
                Orientationx[b, i] = 1.5*(rjk[0]**2) - .5
                a += 1
                b += 1
            a += 2
            
            
            
    if avgT:
        OrientationxAvgT = np.mean(Orientationx,axis=1)
        return OrientationxAvgT
    else:
        return Orientationx
    
    
def Persistence_Length( Positions, N, T, k, chainlength):
    Num_Chains = N/chainlength
    Num_bonds = (chainlength -1)*Num_Chains
    Orientation = np.zeros((T, Num_Chains, Num_Bonds))
    for i in range(T):
        a =0
        for j in range(Num_Chains):
            U0 = Positions[0, a*i] - Positions[0,a*i + 1]
            U0 = np.linalg.norm(U0)
            for k in range(Num_bonds):
                Uk = Position[0, a*i] - Positions[0, a*i  + 1]
                Uk = np.linalg.norm(Uk)
                Orientation[i, j, k] = U0[0]*Uk[0] + U0[1]*Uk[1] + U0[2]*Uk[2]
                a += 1
            a+=2
    
    OrientationAVGT = np.mean(Orientation, axis = 0)
    OrientationAVGChain = np.mean(OrientationAVGT, axis = 0)
    
    return Orientation, OrientationAVGT, OrientationAVGChain
    
                
                
                
    
def Calculate_LP( In_FILE, Restart_In, Restart_out, Num_Steps=10000, Time_Step = 1, Temp= 300, Press = 1):
    """
    Function for calculating the persistence length of Huang's model of P3HT
    """
    os.mkdir('Persistence Length')
    os.chdir('Persistence Length')
    
    Traj_File = 'Persistance_Length.lammpstrj'
    
    File = open(In_FILE, 'w')
    File.write('# Input file for calculating persistence length of P3HT in melt\n\n')
    File.write('read_restart %s\n' % Restart_In)
    File.write('pair_style table linear 187 \n')
    File.write('pair_coeff 1 1 HuangPair2.table P1P1\n')
    File.write('pair_coeff 1 2 HuangPair2.table P1P2\n')
    File.write('pair_coeff 1 3 HuangPair2.table P1P3\n')
    File.write('pair_coeff 2 2 HuangPair2.table P2P2\n')
    File.write('pair_coeff 2 3 HuangPair2.table P2P3\n')
    File.write('pair_coeff 3 3 HuangPair2.table P3P3\n')
    File.write('fix 1 all npt temp %d %d 100 iso %d %d 1000 drag 2\n' % (Temp, Temp, Press, Press))
    File.write('fix 2 all momentum 1 linear 1 1 1\n')
    File.write('thermo_style custom step temp press etotal density\n')
    File.write('thermo 500\n')
    File.write('timestep %d\n' % Time_Step)
    File.write('dump 1 all custom 100 %s id type mol xs ys zs\n' % Traj_File)
    File.write('dump_modify 1 sort id\n')
    File.write('run %d\n' % Num_Steps)
    File.write('write_restart %s\n' % Restart_out)
    
    File.close()
    os.system('mpirun -np 8 lammps <  %s' % In_FILE)
    
    N = NumChains*ChainLength
    T = Num_Steps 
    
    
    Positions = Extract_Trajectory(Traj_File, N, 3)
    
    Corr, CorrT, CorrCh = Persistence_Length( Positions, N, T, 3, ChainLength)
    
    plt.figure()
    plt.plot(CorrCh)
    plt.show()    
    
    return
    
def Strain(In_File, Restart_In, Restart_Out, Stress_Strain_File , Dump = 'none', Strain_Rate = .000001, Num_Steps = 10000, Time_Step = 1, Press = 0, Temp = 423):
    """
    Function for running a uniaxial tensile deformation simulation with LAMMPS
    
    Required Inputs:
    In_File = name of simulation input file
    
    Restart_In = Name of restart file to read in
    
    Restart_Out = Name of restart file to write out
    
    Stress_Strain_File = Name of file to output stress/strain data (Standard Format)
    
    Optional Inputs:
    Strain_Rate = Rate of deformation 
    
    Num_Steps = Number of steps to integrate equations of motion
    
    Time_Step = Time incremant for integrations
    
    Press = externally applied transverse pressure
    
    Temp = Simulation Temperature
    
    """
    
    File = open( In_File, 'w')
    
    File.write('# Input file for running a uniaxial tensile deformation simulation, Filename = %s\n\n' % In_File)
    File.write('read_restart %s\n' % Restart_In)
    File.write('pair_style table linear 187 \n'); File.write('angle_style table linear 200\n')
    File.write('bond_style table linear 200\n')
    File.write('pair_coeff 1 1 HuangPair2.table P1P1\n')
    File.write('pair_coeff 1 2 HuangPair2.table P1P2\n')
    File.write('pair_coeff 1 3 HuangPair2.table P1P3\n')
    File.write('pair_coeff 2 2 HuangPair2.table P2P2\n')
    File.write('pair_coeff 2 3 HuangPair2.table P2P3\n')
    File.write('pair_coeff 3 3 HuangPair2.table P3P3\n');     File.write('bond_coeff 1 HuangBond.table Bond1\n' )
    File.write('bond_coeff 2 HuangBond.table Bond2\n' )
    File.write('bond_coeff 3 HuangBond.table Bond3\n' )
    File.write('angle_coeff 1 HuangAngle.table Angle1\n')
    File.write('angle_coeff 2 HuangAngle.table Angle2\n')
    File.write('angle_coeff 3 HuangAngle.table Angle3\n')
    File.write('angle_coeff 4 HuangAngle.table Angle4\n')
    File.write('variable tmp equal "lx"\n')
    File.write('variable L0 equal ${tmp}\n')
    File.write('variable strainx equal "(lx-v_L0)/v_L0"\n')
    File.write('variable strainy equal "(ly -v_L0)/v_L0"\n')
    File.write('variable strainz equal "(lz -v_L0)/v_L0"\n')
    File.write('variable stress equal "-pxx/10000*1.01325"\n')
    File.write('variable strainrate equal "%f"\n' % Strain_Rate)
    File.write('thermo_style custom step temp lx ly lz pxx density etotal\n')
    File.write('thermo 300\n')
    
    if (Dump != 'none'):
        File.write('dump 1 all custom 200 %s id type mol xs ys zs vx vy vz\n' % Dump)
        
    File.write('reset_timestep 0\n')
    File.write('timestep %d\n' % Time_Step)
    File.write('fix 1 all npt temp %d %d 100 y %d %d 1000 z %d %d 1000 drag 2\n' % ( Temp, Temp, Press, Press, Press, Press))
    File.write('fix 2 all deform 1 x erate %f remap x units box\n' % Strain_Rate)
    File.write('fix def1 all print 1 "${strainx}  ${strainy}  ${strainz}  ${stress}" file %s screen no\n' % Stress_Strain_File)
    File.write('run %d\n' % Num_Steps)
    File.write('write_restart %s\n' % Restart_Out)
    File.write('unfix 1\n')
    File.write('unfix 2\n')
    File.close()
    
    os.system('mpirun -np 8 lammps < %s ' % In_File)
    return
    
    
    
    
    
    
    
    
