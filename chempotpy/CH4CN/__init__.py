'''
CH4CN
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for CH4CN:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  CH4CN_VBMM:           single-state surface of CH4CN, ground state,      P/G

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule:
    The C atom closer to N will be considered as the second last C atom, i.e. the C in CN.

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for CH4CN:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  CH4CN_VBMM:           single-state surface of CH4CN, ground state,
                              availability: potential energy
                              functional form: VBMM
                              corresponding surface in POTLIB: ch4cn-2016.f
                              special emphsize on: CN + CH4 → HCN + CH3 reaction 
                              ref: J. Espinosa-Garcia, C. Range, Y. V. Suleimanov,
                                   "Kinetics study of the CN + CH4 hydrogen abstraction 
                                   reaction based on a new ab initio analytical full-
                                   dimensional potential energy surface",
                                   Phys. Chem. Chem. Phys. 19, 19341-19351 (2017). 

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule:
    The C atom closer to N will be considered as the second last C atom, i.e. the C in CN.


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_nitrogen=0
    n_carbon=0
    natoms=len(geom)

    if natoms!=7:
        print("number of atoms not equal 7 for CH4CN system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='N'):
                 n_nitrogen=n_nitrogen+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
        if n_nitrogen!=1:
            print("number of Nitrogen atoms not equal 1 for CH4CN system")
            sys.exit()
        if n_carbon!=2:
            print("number of Carbon atoms not equal 2 for CH4CN system")
            sys.exit()
        if n_hydrogen!=4:
            print("number of Hydrogen atoms not equal 4 for CH4CN system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for CH4CN_VBMM: input Cartesian should in order of HCHHHCN
    if surface=='CH4CN_VBMM':
        geom_ordered=sorted(geom, key=lambda x: ("H", "C", "N").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "C", "N").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        r1=np.linalg.norm(xyz[4]-xyz[6])
        r2=np.linalg.norm(xyz[5]-xyz[6])
        if r1<r2:
          tmp_row=np.copy(xyz[4])
          xyz[4]=xyz[5]
          xyz[5]=tmp_row
          tmp_row=np.copy(rank_ordered[4])
          rank_ordered[4]=rank_ordered[5]
          rank_ordered[5]=tmp_row
          tmp_row=np.copy(xyz[4])
          xyz[4]=xyz[1]
          xyz[1]=tmp_row
          tmp_row=np.copy(rank_ordered[4])
          rank_ordered[4]=rank_ordered[1]
          rank_ordered[1]=tmp_row
        else:
          tmp_row=np.copy(xyz[4])
          xyz[4]=xyz[1]
          xyz[1]=tmp_row  
          tmp_row=np.copy(rank_ordered[4])
          rank_ordered[4]=rank_ordered[1]
          rank_ordered[1]=tmp_row

    return xyz, rank_ordered


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
