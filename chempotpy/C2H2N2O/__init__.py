'''
C2H2N2O
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for C2H2N2O:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  C2H2N2O_PIPNN:        single-state surface of C2H2N2O, ground state,      P



        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for C2H2N2O:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  C2H2N2O_PIPNN:        single-state surface of C2H2N2O, ground state,
                              availability: potential energy
                              functional form: PIP-NN
                              corresponding surface in POTLIB: C2H2N2O.zip
                              special emphsize on: N2O + C2H2 → Oxadiazole 
                              reaction
                              ref: Y. Liu, J. Li, 
                                   "Quantitative Dynamics of the N2O + C2H2 → 
                                   Oxadiazole Reaction: A Model for 1,3-Dipolar
                                   Cycloadditions", 
                                   ACS Omega 5, 23343-23350 (2020).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_oxygen=0
    n_carbon=0
    n_nitrogen=0
    natoms=len(geom)

    if natoms!=7:
        print("number of atoms not equal 7 for C2H2N2O system")
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
            print("number of Oxygen atoms not equal 1 for C2H2N2O system")
            sys.exit()
        if n_carbon!=2:
            print("number of Carbon atoms not equal 2 for C2H2N2O system")
            sys.exit()
        if n_hydrogen!=2:
            print("number of Hydrogen atoms not equal 2 for C2H2N2O system")
            sys.exit()
        if n_nitrogen!=2:
            print("number of Nitrogen atoms not equal 2 for C2H2N2O system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for C2H2N2O_PIPNN: input Cartesian should in order of HHCCNNO
    if surface=='C2H2N2O_PIPNN':
        geom_ordered=sorted(geom, key=lambda x: ("H", "C", "N", "O").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "C", "N", "O").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

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
