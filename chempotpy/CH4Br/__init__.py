'''
CH4Br
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for CH4Br:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  CH4Br_VBMM_2002:      single-state surface of CH4Br, ground state,      P/G

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for CH4Cl:
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  CH4Br_VBMM_2002:      single-state surface of CH4Br, ground state,
                              availability: potential energy, gradient
                              functional form: VBMM
                              corresponding surface in POTLIB: ch4br.f
                              special emphsize on: CH3 + HBr → CH4 + Br reaction 
                              ref: J. Espinosa-Garcia, 
                                   "Potential energy surface for the CH3 + HBr → 
                                   CH4 + Br hydrogen abstraction reaction: Thermal 
                                   and state-selected rate constants, and kinetic 
                                   isotope effects",
                                   J. Chem. Phys. 117, 2076-2086 (2006).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_bromine=0
    n_carbon=0
    natoms=len(geom)

    if natoms!=6:
        print("number of atoms not equal 6 for CH4Br system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='Br'):
                 n_bromine=n_bromine+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
        if n_bromine!=1:
            print("number of Bromine atoms not equal 1 for CH4Br system")
            sys.exit()
        if n_carbon!=1:
            print("number of Carbon atoms not equal 1 for CH4Br system")
            sys.exit()
        if n_hydrogen!=4:
            print("number of Hydrogen atoms not equal 4 for CH4Br system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for CH4Br_VBMM_2002: input Cartesian should in order of HCHHHBr
    if surface=='CH4Br_VBMM_2002':
        geom_ordered=sorted(geom, key=lambda x: ("C", "H", "Br").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("C", "H", "Br").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        tmp_row=np.copy(xyz[1])
        xyz[1]=xyz[0]
        xyz[0]=tmp_row
        tmp_row=np.copy(rank_ordered[1])
        rank_ordered[1]=rank_ordered[0]
        rank_ordered[0]=tmp_row

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
