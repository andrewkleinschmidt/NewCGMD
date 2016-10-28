import numpy as np
import math 
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import quad
import sys
import Huangparameters as Huang

def Avg_Val_Flory_Schulz(alpha, x):
    return alpha**2 * x**2 * (1-alpha)**(x-1)

def Flory_Schulz(alpha,x):
    return alpha**2 * x * (1-alpha)**(x-1)

def Gaussian(sigma,x,Mn):
    return math.sqrt(2 * sigma **2 * math.pi) ** -1 * math.exp(-1 * ((x - Mn)**2) / (2 * sigma **2))

def find_polydispersity(sigma, mean):
    Mn = 0
    Mw = 0
    for i in range(0,10*mean):
        Mn += Gaussian(sigma,i,mean)*i
        Mw += Gaussian(sigma,i,mean)*i**2
    Mw /= Mn
    return (Mw/Mn)

"""
    For flory-schulz
    def find_polydispersity(alpha, avg):
    Mn = 0
    Mw = 0
    for i in range(0,10*avg):
        Mn += Flory_Schulz(alpha,i)*i
        Mw += Flory_Schulz(alpha,i)*i**2
        print Mn
        print Mw
    print Mn
    print Mw
    Mw /= Mn
    return (Mw/Mn)"""

def optimize_sigma(deg_polymerization,dispersity):
    print "I ran"
    lowbound = 0.
    highbound = float(deg_polymerization) * 10.
    sigma = float((highbound - lowbound)/2)
    est_dispersity = find_polydispersity(sigma,deg_polymerization)
    print est_dispersity
    print dispersity
    while abs(est_dispersity-dispersity)/dispersity > .001:
        if est_dispersity < dispersity:
            lowbound, sigma = sigma, (sigma + highbound)/2
        else:
            highbound, sigma = sigma, (sigma + lowbound)/2
        est_dispersity = find_polydispersity(sigma,deg_polymerization)
        print sigma
    return sigma

def make_bins(sigma,deg_polymerization,num_polymer):
    """ This function creates a binned distribution of polymers fitting a Gaussian curve.  It first generates the probabilities over a range of degrees of polymerization,
        then sums them in a bin size (currently given as 50) to find the average total probability of finding a polymer around that size.  It only adds bins of at least a certain
        total probability (here at least 0.1%), and reweights at the end to give a total probability of 1.
        dx: small step across degrees of polymerization, x: all x values to sum over, PDF: the correspoonsing Gaussian probability data to x, bin: temp variable to add all terms in 
        a bin while summing, bin_list: list of bin probabilities over a certain threshold, index_list: list of corresponding degrees of polymerization for bin probabilities."""
    dx = .01
    x = np.arange(deg_polymerization*10, step = dx)
    PDF = []
    for i in range(len(x)):
        PDF.append(Gaussian(sigma,x[i],deg_polymerization))
    bin = 0
    bin_list = []
    index_list = []
    for i in range (1,len(PDF)):
        bin += dx*PDF[i]
        if (dx * i) % 50 == 0:
            if bin > .01:
                print dx*i
                bin_list.append(bin)
                index_list.append(int((dx*i*2-50)/2))
            bin = 0
    total = 0
    for b in bin_list:
        total += b
    bin_list[:] = [int(b * 1/total * num_polymer) for b in bin_list]
    """Mn = sum([a*b for a,b in zip(index_list,bin_list)])
    Mw = sum([a**2*b for a,b in zip(index_list,bin_list)])
    Mw /= Mn <- can be used to check polydispersity of sample correct"""
    return bin_list, index_list

def gen_CDF_for_size_selection(alpha, Mw):
    dx = .01
    x = np.arange(0, Mw*50)
    PDF = np.fromfunction(Flory_Schulz(alpha,x), (len(x),1))
    CDF = np.full_like(x, 0.0)
    print PDF
    for i in range (1,len(PDF)):
        CDF[i] = CDF[i-1] + dx*PDF[i]
    return CDF

def select_polymers(CDF, num_polymer):
    polymers = []
    random.seed()
    for i in range (0,num_polymer):
        rand = random.random()
        j = 0
        while CDF[j]<rand:
            j+=1
        polymers.append(int(j))
    return polymers

def Gen_PDF_CDF_OPLS( V, Beta ):
	"""
	This function takes in a numpy array V that containes the energetic coefficients for the OPLS style dihedral potential of the form:
			U = (1/2)V1(1+cos(phi)) + (1/2)V2(1-cos(2phi)) + (1/2)V3(1+cos(3phi)) + ....
	It then uses Boltzmann statistics along with the inverse temperature Beta to generate a PDF and CDF of the dihedral angle
	
	The output is two numpy arrays that represent the PDF and CDF associated with this potential energy function
	"""
	dx = 0.0001
	x = np.arange(0, 6.28, dx) # Generate discretized values for x (phi)
	U = np.full_like(x, 0.0) # Initialize potential energy array
	PDF = np.full_like(x, 0.0) # Initialize PDF array
	CDF_NN = np.full_like(x, 0.0) # Initialize non-normalized CDF array
	CDF = np.full_like(x, 0.0) # Initialize normalized CDF array
	norm = 0
	L = len(x.tolist()) 
	U = 0.5*(V[0]*(1 + np.cos(x)) + V[1]*(1 - np.cos(2*x)) + V[2]*(1 + np.cos(3*x)) + V[3]*(1 - np.cos(4*x)))
	PDF = np.exp(-U*Beta)
	
	for i in range(L-1):
		CDF_NN[i+1] = CDF_NN[i] + PDF[i]*dx
	
	for i in range(L):
		PDF[i] = PDF[i]/CDF_NN[-1]
		norm += PDF[i]*dx
	
	for i in range(L-1):
		CDF[i+1] = CDF[i] + PDF[i]*dx
		
	return PDF, CDF
 
def Gen_PDF_CDF_Multi( V, Beta ):
	"""
	This function takes in a numpy array V that containes the energetic coefficients for the Multi style dihedral potential of the form
	It then uses Boltzmann statistics along with the inverse temperature Beta to generate a PDF and CDF of the dihedral angle
	
	The output is two numpy arrays that represent the PDF and CDF associated with this potential energy function
    
    Changed to only return CDF
	"""
	dx = 0.0001
	x = np.arange(0, 6.28, dx) # Generate discretized values for x (phi)
	U = np.full_like(x, 0.0) # Initialize potential energy array
	PDF = np.full_like(x, 0.0) # Initialize PDF array
	CDF_NN = np.full_like(x, 0.0) # Initialize non-normalized CDF array
	CDF = np.full_like(x, 0.0) # Initialize normalized CDF array
	norm = 0
	L = len(x.tolist()) 
	U = V[0] + V[1]*np.cos(x) + V[2]*np.cos(x)**2 + V[3]*np.cos(x)**3 + V[4]*np.cos(x)**4
	PDF = np.exp(-U*Beta)
	
	for i in range(L-1):
		CDF_NN[i+1] = CDF_NN[i] + PDF[i]*dx
	
	for i in range(L):
		PDF[i] = PDF[i]/CDF_NN[-1]
		norm += PDF[i]*dx
	
	for i in range(L-1):
		CDF[i+1] = CDF[i] + PDF[i]*dx
		
	return CDF
	
def Plot_Dihedral( V, style):
    """ 
    This function plots the dihedral energy function from 0 to 2*pi
    inputs:
        V is a numpy array containing the parameters for the functional form of the dihedral potential
        style can either be 'OPLS' or 'MULTI'
    returns nothing
    """
    dx = 0.0001
    x = np.arange(0, 6.28, dx)
    U = np.full_like(x, 0.0)    
    
    if style == 'OPLS':
        U = 0.5*(V[0]*(1 + np.cos(x)) + V[1]*(1 - np.cos(2*x)) + V[2]*(1 + np.cos(3*x)) + V[3]*(1 - np.cos(4*x)))
        
    elif style == 'MULTI':
        print 'Under Construction! Sorry'
        
    
    plt.figure()
    plt.plot(x, U)
    plt.title('Dihedral Potential')
    plt.xlabel('Dihedral Angle (rad)')
    plt.ylabel('Energy (kcal/mol)')
    plt.xlim((0,6.28))
    plt.show()
    return
    
def Compare_Bond( Bond ):
    dx = 0.0001
    N = len(Bond) - 1
    r = np.arange(Bond[0]-.5, Bond[0]+.5, dx)
    U1 = np.full_like(r, 0.0)
    U2 = np.full_like(r, 0.0)
    
    for i in range(1,N+1):
        U1 += Bond[i]*(r - Bond[0])**(i+1)
        print i+1
        
    U2 = Bond[1]*(r - Bond[0])**2
    
    plt.figure()
    plt.plot(r, U1, label = 'Full Potential')
    plt.plot(r, U2, label = 'Harmonic Approximation')
    plt.xlim((Bond[0]-.5,Bond[0]+.5))
    plt.ylim((U1.min(), U1.max()))
    plt.title('P2P3', fontsize = 30)
    plt.xlabel('Bond Length (Angstrom)', fontsize = 20)
    plt.ylabel('Potential Energy (Kcal/mol)', fontsize = 20)
    plt.legend()
    plt.show()
    return
    
def Compare_Angle( Angle):
    dx = 0.0001
    M = [0, 2 , 3, 4 , 5 , 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
    N = len(Angle) - 1
    Th0 = Angle[0]*(3.1415/180.)
    dTh = .5
    Th = np.arange(Th0-dTh, Th0 + dTh, dx)
    U1 = np.full_like(Th, 0.0)
    U2 = np.full_like(Th, 0.0)
    
    for i in range(1, N+1):
        U1 += Angle[i]*(Th - Th0)**M[i-1]
    
    U2 = Angle[2]*(Th - Th0)**2 + Angle[3]*(Th - Th0)**3 + Angle[4]*(Th-Th0)**4
    plt.figure()
    plt.plot(Th, U1, label = 'Full Potential')
    plt.plot(Th, U2, label = 'Class 2 Approximation')
    plt.ylim((U2.min(), U1.max() ))
    plt.xlim((Th0-dTh, Th0 + dTh))
    plt.title('P2P1P1')
    plt.xlabel('Angle (Radians)', fontsize = 20)
    plt.ylabel('Potential Energy (Kcal/mol)', fontsize = 20)
    plt.legend()
    plt.show()

    return  

def Compare_Dihedral(Dihedral):
    dx = 0.0001
    N = len(Dihedral)
    Dih = np.arange(0, 2*3.1415, dx)
    U1 = np.full_like(Dih, 0.0)
    U2 = np.full_like(Dih, 0.0)
    
    for i in range(N):
        U1 += Dihedral[i]*np.cos(Dih)**i
    
    for i in range(4):
        U2 += Dihedral[i]*np.cos(Dih)**i
    
    plt.figure()
    plt.plot(Dih, U1, label = 'Full Potential')
    plt.plot(Dih, U2, label = 'Truncated Approximation')
    plt.ylim((U2.min(), U2.max()))
    plt.xlim((0.0, 2*3.1415))
    plt.title('P2P1P1P2', fontsize = 30)
    plt.xlabel('Dihedral Angle (radians)', fontsize = 20)
    plt.ylabel('Energy (Kcal/mol)', fontsize = 20)
    plt.legend()
    plt.show()
    return
    
            
        
    

def Gen_Random_Dih(CDF):
    """
    Function for generating a random Dihedral angle by inverting the CDF
    input: Numpy array CDF
    Output: Random Float [0,2*pi]
    """
    dx = 0.0001
    L = len(CDF.tolist())
    h = random.random()
    for i in range(L-1):
        if h > CDF[i]:
            return i * dx
    return L - 1
 
def Normalize(A):
    """
    Function for normalizing a Numpy 3-array
    """
    norm = 0
    for i in range(len(A.tolist())):
        norm += A[i]**2
    B = A/(math.sqrt(norm))
    return B
	
def dotproduct(A,B):
	# Function to compute dotproduct
	C = A[0]*B[0] + A[1]*B[1] + A[2]*B[2]
	return C
	
def crossproduct(A,B):
	# Function to compute cross product
	C= np.zeros(3,float)
	C[0] = A[1]*B[2] - A[2]*B[1]
	C[1] = A[2]*B[0] - A[0]*B[2]
	C[2] = A[0]*B[1] - A[1]*B[0]
	return C
 
 

def Apply_PBC(Position, Box_Length, substrate = False):
    """ Apply Periodic Boundary Conditions to a set of particle positions
        
        input: 
               Position: Numpy array [N,3] Containing position coordinates of particles
               
               Box_Length: Length of Simulation box (float)
               
        Output: 
               PositionPBC: Numpy array [N,3] Containing remapped coordinates of particles
        
        If there is a LJ wall representing the substrate, do not create periodic boundaries for z
    """
    N = np.size(Position,0)
    if substrate:
        M = np.size(Position,1) - 1
    else:
        M = np.size(Position,1)
    PositionPBC = np.zeros((N,3), float)
    for i in range(N):
        for j in range(M):
            PositionPBC[i,j] = Position[i,j] - Box_Length*math.floor(Position[i,j]/Box_Length)
        if substrate:
            PositionPBC[i,2] = Position[i,2]
    return PositionPBC

def Gen_rand_orthonormal_Basis():
    """
    This function generates a random orthonormal basis in the form of a 3x3 numpy array
    """
    xrsq = random.random()
    yrsq= (1 - xrsq)*random.random()
    zrsq = 1 - xrsq -yrsq
    yr = (random.randint(0,1)*2 - 1)*math.sqrt(yrsq)
    zr = (random.randint(0,1)*2 - 1)*math.sqrt(zrsq)
    xr = (random.randint(0,1)*2 - 1)*math.sqrt(xrsq)
    Basis = np.zeros((3,3))
    Basis[0] = [xr, yr, zr]
    x2r = random.random()
    y2r = -xr*x2r*random.random()
    z2r = -(xr*x2r + yr*y2r)/zr
    Basis[1] = [x2r,y2r, z2r]
    Basis[2] = crossproduct(Basis[0], Basis[1])
    return Basis
    
def Gen_Linear_Polymer(ChainLength, r0, th0, CDF, Box_Length):
    Bases = np.zeros([ChainLength,3,3], float)
    Position = np.zeros([ChainLength,3], float)
    Position[0] = [ Box_Length*random.random(), Box_Length*random.random(), Box_Length*random.random()]
    Bases[0] = Gen_rand_orthonormal_Basis()
    
    
    for i in range(ChainLength-1):
        Phi = Gen_Random_Dih(CDF)
        Position[i+1] = Position[i]
        Position[i+1] += r0*(math.cos(th0+3.14))*Bases[i,0]
        Position[i+1] += r0*(math.sin(th0+3.14))*(math.cos(Phi))*Bases[i,1]
        Position[i+1] += r0*(math.sin(th0))*(math.sin(Phi))*Bases[i,2]
        Bases[i+1, 0] = Normalize(Position[i+1] - Position[i])
        Bases[i+1, 1] = Normalize(crossproduct(Bases[i+1,0], Bases[i,0]))
        Bases[i+1, 2] = crossproduct(Bases[i+1,0], Bases[i+1,1])
        
    Position = Apply_PBC(Position, Box_Length)
    
    return Position, Bases
    
def generatepolymer( chainlength, BondLengths , Angles, CUM2, Box_Length, zbox_length, substrate = False):
    Bases = np.zeros([chainlength,3,3],float)
    Position_M = np.zeros([chainlength,3],float)
    Position_S1 = np.zeros([chainlength,3],float)
    Position_S2 = np.zeros([chainlength,3],float)
    Position_M[0] = [ Box_Length*random.random(), Box_Length*random.random(), zbox_length*random.random()]
    #Bases[0] = [[1,0,0],[0,1,0],[0,0,1]]
    Bases[0] = Gen_rand_orthonormal_Basis()
    Position_S1[0] = Position_M[0] + BondLengths['P1P2']*Bases[0,2]
    Position_S2[0] = Position_S1[0] + BondLengths['P2P3']*Bases[0,2]
    for i in range(chainlength-1):
        R = BondLengths['P1P1']
        Theta = 0.3 - 0.6 * random.random()
        Phi = Gen_Random_Dih(CUM2)
        Position_M[i+1] = Position_M[i]
        Position_M[i+1] += R*(math.cos(Theta))*Bases[i,0] 
        Position_M[i+1] += R*(math.cos(Phi))*(math.sin(Theta))*Bases[i,1]
        Position_M[i+1] += R*(math.sin(Theta))*(math.sin(Phi))*Bases[i,2]
        Bases[i+1,0] = Normalize(Position_M[i+1] - Position_M[i])
        Bases[i+1,1] = Normalize(crossproduct(Bases[i+1,0],Bases[i,0]))
        Bases[i+1,2] = crossproduct(Bases[i+1,0],Bases[i+1,1])
        Position_S1[i+1] = Position_M[i+1] + BondLengths['P1P2']*Bases[i+1,2]
        Position_S2[i+1] = Position_S1[i+1] + BondLengths['P2P3']*Bases[i+1,2]
        #print dotproduct(Bases[i+1,0], Bases[i+1,1])
            #print dotproduct(Bases[i+1,0], Bases[i+1,2])
		#print dotproduct(Bases[i+1,1], Bases[i+1,2])
    return Position_M, Position_S1, Position_S2, Bases

    
    
def Quiver_Plot(Position, Bases, Box_Length):
    """
	Function that takes an array of positions and their associated orthonormal basis sets to 
	generate a 3D vector plot of the polymer chains along with their respective orientations
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection= '3d')
    ax.quiver(Position[:,0], Position[:,1], Position[:,2], Bases[:,0,0], Bases[:,0,1], Bases[:,0,2], length = 3, color='r')
    ax.quiver(Position[:,0], Position[:,1], Position[:,2], Bases[:,1,0], Bases[:,1,1], Bases[:,1,2], length = 3, color='g')
    ax.quiver(Position[:,0], Position[:,1], Position[:,2], Bases[:,2,0],Bases[:,2,1], Bases[:,2,2], length = 3)
    ax.scatter(Position[:,0], Position[:,1], Position[:,2],'k')
    #ax.grid(on=False)
    plt.title('Non-Overlapping Random Configuration, PBC\'s')
    plt.xlim((0,Box_Length))
    plt.ylim((0,Box_Length))
    plt.show()
    return



def Check_For_Overlap(Polymer_Coords, New_Poly, radius):
    """ 
    Function to check for overlap between two polymer chains
    input: Polymer_Coords = Existing polymer chains
            New_Poly = Potential New Polymer Chain
            radius = radius... DUH
    Output: True indicates no overlap
            False indicates overlap
    """
    Rij = np.zeros(3)
    radiusSQ = radius*radius
    
    for polymer in Polymer_Coords:
        midpoint = np.size(polymer,0)/2
        new_mid = np.size(New_Poly,0)/2
        R00 = polymer[0] - New_Poly[0]
        R00 = math.sqrt(dotproduct(R00,R00))
        Rmidmid = polymer[midpoint] - New_Poly[new_mid]
        Rmidmid = math.sqrt(dotproduct(Rmidmid,Rmidmid))
        Rlastlast = polymer[-1] - New_Poly[-1]
        Rlastlast = math.sqrt(dotproduct(Rlastlast,Rlastlast))
        if R00 >1000 and Rmidmid > 1000 and Rlastlast > 1000:
            continue
        for i in range(np.size(polymer,0)):
            for j in range(np.size(New_Poly,0)):
                Rij = polymer[i] - New_Poly[j]
                RijSQ = dotproduct(Rij,Rij)
                if (RijSQ < radiusSQ):
                    print 'Overlap'
                    return False
    return True     
    
def Gen_Many_Lin_Polymers( ChainLength, NumChains, r0, th0, CDF, Sigma, Box_Length):
    """ 
    Function to generate many non-overlapping polymer chains
    """
    Position, Bases = Gen_Linear_Polymer(ChainLength, r0, th0, CDF,  Box_Length)
    NChains = 1
    while (NChains < NumChains):
        PositionNew, BasesNew = Gen_Linear_Polymer(ChainLength, r0, th0, CDF, Box_Length)
        Bool = Check_For_Overlap(Position, PositionNew, Sigma)
        if (Bool):
            Position = Apply_PBC( np.concatenate((Position, PositionNew), axis=0), Box_Length)
            Bases = np.concatenate((Bases, BasesNew), axis=0)
            NChains += 1
            print "Polymer Deposited %d " % NChains
    return Position, Bases
    


def Gen_Many_Polymers( AvgChainLength, NumChains, BondLengths, Angles, SigmaM_M, Box_Length, zbox_length,Dispersity, substrate = False):
    """
    Function to generate many non-overlapping polymer chains and return modified z box size
    """
    CUM2 = Gen_PDF_CDF_Multi(Huang.P1P1P1P1, Huang.Beta)
    num_polys, lengths = make_bins(optimize_sigma(AvgChainLength,Dispersity),AvgChainLength,NumChains)
    zhi = 0
    zlo = 0
    num = 0
    Position_M, Position_S1, Position_S2, Bases = generatepolymer( lengths[0], BondLengths , Angles, CUM2, Box_Length, zbox_length, substrate)
    Polymer_Coords_M = []
    Polymer_Coords_S1 = []
    Polymer_Coords_S2 = []
    Polymer_Coords_M.append(Apply_PBC(Position_M,Box_Length,substrate))
    Polymer_Coords_S1.append(Apply_PBC(Position_S1,Box_Length,substrate))
    Polymer_Coords_S2.append(Apply_PBC(Position_S2,Box_Length,substrate))
    for i in range(np.size(Position_M,0)):
        if zhi < Position_M[i][2]:
            zhi = Position_M[i][2]
        if zhi < Position_S1[i][2]:
            zhi = Position_S1[i][2]
        if zhi < Position_S2[i][2]:
            zhi = Position_S2[i][2]
        if zlo > Position_M[i][2]:
            zlo = Position_M[i][2]
        if zlo > Position_S1[i][2]:
            zlo = Position_S1[i][2]
        if zlo > Position_S2[i][2]:
            zlo = Position_S2[i][2]
    for i in range(len(lengths)):
        for q in range(num_polys[i]):
                num += 1
            #Bool = False
            #while Bool == False:
                PositionNew_M, PositionNew_S1, PositionNew_S2, BasesNew = generatepolymer( lengths[i], BondLengths , Angles, CUM2, Box_Length,zbox_length, substrate)
                #Bool = Check_For_Overlap(Polymer_Coords_M, PositionNew_M, SigmaM_M)
                #if Bool:
                Polymer_Coords_M.append(PositionNew_M)
                Polymer_Coords_S1.append(PositionNew_S1)
                Polymer_Coords_S2.append(PositionNew_S2)
                Bases = np.concatenate((Bases, BasesNew), axis=0)
                print "Polymer Deposited"
                print num
                for h in range(np.size(PositionNew_M,0)):
                    if zhi < PositionNew_M[h][2]:
                        zhi = PositionNew_M[h][2]
                    if zhi < PositionNew_S1[h][2]:
                        zhi = PositionNew_S1[h][2]
                    if zhi < PositionNew_S2[h][2]:
                        zhi = PositionNew_S2[h][2]
                    if zlo > PositionNew_M[h][2]:
                        zlo = PositionNew_M[h][2]
                    if zlo > PositionNew_S1[h][2]:
                        zlo = PositionNew_S1[h][2]
                    if zlo > PositionNew_S2[h][2]:
                        zlo = PositionNew_S2[h][2]
    print zhi
    print zlo
    for i in range(np.size(Polymer_Coords_M,0)):
        for j in range(np.size(Polymer_Coords_M[i],0)):
            Polymer_Coords_M[i][j][2] += -1 * zlo + .001
            Polymer_Coords_S1[i][j][2] += -1 * zlo + .001
            Polymer_Coords_S2[i][j][2] += -1 * zlo + .001
    zBoxLengthNew = zhi - zlo + 1

    return Polymer_Coords_M, Polymer_Coords_S1, Polymer_Coords_S2, Bases, zBoxLengthNew
    
	
def Deposit_PCBM( Position_M, SigmaM_P, SigmaP_P, NumPCBM, Box_Length):
    """
    Function to deposit PCBM Molecules randomly so that they dont overlap with the polymer molecules
    """
    
    NS = 0
    Pos = np.zeros(3,float)
    Position_C60 = np.zeros((NumPCBM,3), float)
    NumM = np.size(Position_M,0)
    SigmaM_PSQ = SigmaM_P*SigmaM_P
    SigmaP_PSQ = SigmaP_P*SigmaP_P
    while (NS<NumPCBM):
         Pos = [Box_Length*random.random() - 20, Box_Length*random.random() - 20,Box_Length*random.random() - 20 ]
         deposit = True
         for i in range(NumM):
             RijMP = Pos - Position_M[i]
             RijSQMP = dotproduct(RijMP,RijMP)
             if (RijSQMP < SigmaM_PSQ):
                 deposit = False
 		 print "Fullerene-P3HT Overlap"
         if (deposit and NS>0):
            for i in range(NS):
                RijPP = Pos - Position_C60[i]
                RijPPSQ = dotproduct(RijPP,RijPP)
                if (RijPPSQ < SigmaP_PSQ):
                    deposit = False
                    print "Fullerene-Fullerene Overlap"
        
         if (deposit):
            Position_C60[NS] = Pos
            NS += 1
            print "PCBM Deposited"    
    return Position_C60


            
