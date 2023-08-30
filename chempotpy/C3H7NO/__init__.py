'''
C3H7NO
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for C3H7NO:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  C3H7NO_PIP:           single-state surface of C3H7NO, ground state,      P/G

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule:
    The structure is: CH3-NH-CO-CH3. Because this surfaces only performs partially 
    permutationally invariant, therefore, we compute NH, CO, CCO, CH distances to 
    figure out the ordering. 


        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for C3H7NO:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  C3H7NO_PIP:           single-state surface of C3H7NO, ground state,
                              availability: potential energy
                              functional form: PIP
                              corresponding surface in POTLIB: C3H7NO.zip
                              ref: A. Nandi, C. Qu, and J. M. Bowman,
                                   "Full and fragmented permutationally invariant
                                   polynomial potential energy surfaces for trans
                                   and cis N-methyl acetamide and isomerization 
                                   saddle points",
                                   J. Chem. Phys. 151, 084306 (2019).

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule:
    The structure is: CH3-NH-CO-CH3. Because this surfaces only performs partially 
    permutationally invariant, therefore, we compute NH, CO, CCO, CH distances to 
    figure out the ordering. 

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_oxygen=0
    n_carbon=0
    n_nitrogen=0
    natoms=len(geom)

    if natoms!=12:
        print("number of atoms not equal 12 for C3H7NO system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                n_oxygen=n_oxygen+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='N'):
                 n_nitrogen=n_nitrogen+1
        if n_oxygen!=1:
            print("number of Oxygen atoms not equal 1 for C3H7NO system")
            sys.exit()
        if n_carbon!=3:
            print("number of Carbon atoms not equal 3 for C3H7NO system")
            sys.exit()
        if n_hydrogen!=7:
            print("number of Hydrogen atoms not equal 7 for C3H7NO system")
            sys.exit()
        if n_nitrogen!=1:
            print("number of Nitrogen atoms not equal 1 for C3H7NO system")
            sys.exit()

    xyz=np.zeros((natoms,3))
   
    rank=np.copy(geom)
    ordered_rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for C3H7NO_PIP: input Cartesian should in order of HHHCNHOCCHHH
    if surface=='C3H7NO_PIP':
        geom_ordered=sorted(geom, key=lambda x: ("H", "C", "N", "O").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "C", "N", "O").index(x[0]))
        #geom_ordered[9]=tmp_row
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        array_H=[0, 1, 2, 3, 4, 5, 6]
        array_C=[7, 8, 9]
        #first, figrue out which H is associated to N
        rNH=np.zeros((7))
        for i in range(7):
          rNH[i]=np.linalg.norm(xyz[i]-xyz[10])
        minNH_idx=np.argmin(rNH)
        row_H_N=xyz[minNH_idx].copy()
        array_H.remove(minNH_idx)
        #second, figure out which C is associated to O 
        rCO=np.zeros((3))
        for i in range(3):
            rCO[i]=np.linalg.norm(xyz[7+i]-xyz[11])
        minCO_idx=np.argmin(rCO)+7
        array_C.remove(minCO_idx)
        #third, figure out which C is associated to CO
        rCC=np.zeros((2))
        for i in range(len(array_C)):
            rCC[i]=np.linalg.norm(xyz[array_C[i]]-xyz[minCO_idx])
        minCC_idx=array_C[np.argmin(rCC)]
        array_C.remove(minCC_idx)
        #now do CH for the C that is associated to CO
        rCH=np.zeros((6))
        for i in range(len(array_H)):
            rCH[i]=np.linalg.norm(xyz[array_H[i]]-xyz[minCC_idx])
        minCH_1, minCH_2, minCH_3, minCH_1_idx, minCH_2_idx, minCH_3_idx=find_three_smallest_numbers(rCH)
        CH1_idx=array_H[minCH_1_idx]
        CH2_idx=array_H[minCH_2_idx]
        CH3_idx=array_H[minCH_3_idx]
        #array_H.remove(minCH_1_idx)
        #array_H.remove(minCH_2_idx)
        #array_H.remove(minCH_3_idx)
        array_H.remove(CH1_idx)
        array_H.remove(CH2_idx)
        array_H.remove(CH3_idx)
        ordered_xyz=np.zeros((natoms,3))
        k=0
        for i in array_H:
            ordered_xyz[k]=np.copy(xyz[i])
            ordered_rank[k]=np.copy(rank_ordered[i])
            k=k+1
        for i in array_C:
            ordered_xyz[k]=np.copy(xyz[i])
            ordered_rank[k]=np.copy(rank_ordered[i])
        ordered_xyz[4]=np.copy(xyz[10])
        ordered_xyz[5]=np.copy(xyz[minNH_idx])
        ordered_xyz[6]=np.copy(xyz[11])
        ordered_xyz[7]=np.copy(xyz[minCO_idx])
        ordered_xyz[8]=np.copy(xyz[minCC_idx])
        ordered_xyz[9]=np.copy(xyz[CH1_idx])
        ordered_xyz[10]=np.copy(xyz[CH2_idx])
        ordered_xyz[11]=np.copy(xyz[CH3_idx])
        xyz=np.copy(ordered_xyz)

        ordered_rank[4]=np.copy(rank_ordered[10])
        ordered_rank[5]=np.copy(rank_ordered[minNH_idx])
        ordered_rank[6]=np.copy(rank_ordered[11])
        ordered_rank[7]=np.copy(rank_ordered[minCO_idx])
        ordered_rank[8]=np.copy(rank_ordered[minCC_idx])
        ordered_rank[9]=np.copy(rank_ordered[CH1_idx])
        ordered_rank[10]=np.copy(rank_ordered[CH2_idx])
        ordered_rank[11]=np.copy(rank_ordered[CH3_idx])
        rank_ordered=np.copy(ordered_rank)


    return xyz, rank_ordered

def find_three_smallest_numbers(arr):
    min1 = arr[0]
    min2 = float('inf')
    min3 = float('inf')
    min1_idx = 0
    min2_idx = -1
    min3_idx = -1
    
    for i in range(1, len(arr)):
        if arr[i] < min1:
            min3 = min2
            min3_idx = min2_idx
            min2 = min1
            min2_idx = min1_idx
            min1 = arr[i]
            min1_idx = i
        elif arr[i] < min2 and arr[i] != min1:
            min3 = min2
            min3_idx = min2_idx
            min2 = arr[i]
            min2_idx = i
        elif arr[i] < min3 and arr[i] != min1 and arr[i] != min2:
            min3 = arr[i]
            min3_idx = i
            
    return min1, min2, min3, min1_idx, min2_idx, min3_idx

def reverse_order(v, g, d, geom, rank_ordered):

    natoms=len(geom)
    nstates=len(v)
    v_ro=np.zeros((nstates))
    g_ro=np.zeros((nstates,natoms,3))
    d_ro=np.zeros((nstates,nstates,natoms,3))

    v_ro = v
    for istate in range(nstates):
        for iatom in range(natoms):
            for idir in range(3):
                g_ro[istate][int(rank_ordered[iatom][1])][idir]=g[istate][iatom][idir]

    for istate in range(nstates):
        for jstate in range(nstates):
            for iatom in range(natoms):
                for idir in range(3):
                    d_ro[istate][jstate][int(rank_ordered[iatom][1])][idir]=d[istate][jstate][iatom][idir]


    return v_ro, g_ro, d_ro
