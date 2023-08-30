'''
CH3NH2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for CH3NH2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Multi-State Surfaces:
    1.  CH3NH2_APRP:       multi-state surface of CH3NH2, 1-2 state,      P/G/D
    2.  CH3NH2_APRP_DPEM:  multi-state surface of CH3NH2, 1-2 state,      U/UG

    ==VERY IMPORTANT NOTICE==
    The atomic ordering specified in geom is required for CH3NH2_APRP
    Internal numbering of atoms:
        H5         H1 
          \       / 
        H6--C4---N3 
          /       \  
        H7         H2 

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for CH3NH2:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 
    APRP: anchor points reactive potential

    Multi-State Surfaces:
    1.  CH3NH2_APRP:       multi-state surface of CH3NH2, 1-3 state,
                           availability: potential energy, gradient, nonadiabatic coupling vector 
                           functional form: APRP + PIP
                           corresponding surface in POTLIB: CH3NH2.f
                           ref: K. A. Parker, and D. G. Truhlar, 
                                "Semiglobal diabatic potential energy matrix for 
                                the N–H photodissociation of methylamine",
                                J. Chem. Phys. 152, 244309 (2020).

    VERY IMPORTANT NOTICE:
    The atomic ordering specified in geom is required for CH3NH2_APRP
        H5         H1  
          \       / 
        H6--C4---N3 
          /       \ 
        H7         H2 

        """)

    return


def check(system, surface, geom):

    n_carbon=0
    n_hydrogen=0
    n_nitrogen=0
    natoms=len(geom)

    if natoms!=7:
        print("number of atoms not equal 7 for CH3NH2 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='N'):
                n_nitrogen=n_nitrogen+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
        if n_nitrogen!=1:
            print("number of Nitrogen atoms not equal 1 for CH3NH2 system")
            sys.exit()
        if n_hydrogen!=5:
            print("number of Hydrogen atoms not equal 6 for CH3NH2 system")
            sys.exit()
        if n_carbon!=1:
            print("number of Carbon atoms not equal 6 for CH3NH2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for CH3NH2_APRP: input Cartesian should in order of HHNCHHH
    if surface=='CH3NH2_APRP' or surface=='CH3NH2_APRP_DPEM':
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom[iatom][idir+1]

    return xyz, rank

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

def reverse_order_dpem(u, ug, geom, rank_ordered):

    natoms=len(geom)
    nstates=len(u)
    u_ro=np.zeros((nstates,nstates))
    ug_ro=np.zeros((nstates,nstates,natoms,3))

    u_ro = u

    for istate in range(nstates):
        for jstate in range(nstates):
            for iatom in range(natoms):
                for idir in range(3):
                    ug_ro[istate][jstate][int(rank_ordered[iatom][1])][idir]=ug[istate][jstate][iatom][idir]


    return u_ro, ug_ro
