'''
C7H8S
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for C7H8S:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Multi-State Surfaces:
    1.  C7H8S_APRP:           multi-state surface of C7H8S, 1-3 state,      P/G/D
    2.  C7H8S_III_APRP:       multi-state surface of C7H8S, 1-3 state,      P/G/D
    3.  C7H8S_S22_APRP:       multi-state surface of C7H8S, 1-3 state,      P/G/D
    4.  C7H8S_APRP_DPEM:      multi-state surface of C7H8S, 1-3 state,      U/UG
    5.  C7H8S_III_APRP_DPEM:  multi-state surface of C7H8S, 1-3 state,      U/UG
    6.  C7H8S_S22_APRP_DPEM:  multi-state surface of C7H8S, 1-3 state,      U/UG


    ==VERY IMPORTANT NOTICE==
    The atomic ordering specified in geom is required for all APRP surfaces for C7H8S
    Internal numbering of atoms:
         H10      H9   H14,15
          \      /      \ 
           C5---C4       C13---H16
          /      \      /
    H11--C6       C3---S12 
          \      /      
           C1---C2      
          /      \ 
         H7       H8

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for C7H8S:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 
    APRP: anchor points reactive potential

    Multi-State Surfaces:
    1.  C7H8S_APRP:           multi-state surface of C7H8S, 1-3 state,
                              availability: potential energy, gradient, nonadiabatic coupling vector 
                              functional form: APRP + PIP
                              corresponding surface in POTLIB: phsh3_aprp.f90
                              pecial emphsize on: C6H5SCH3 → C6H5S + CH3
                              ref: L. Li, and D. G. Truhlar, 
                                   "Full-dimensional ground- and excited-state potential energy
                                   surfaces and state couplings for photodissociation of thioanisole",
                                   J. Chem. Phys. 146, 064301 (2017).
    ============================================================================================
    2.  C7H8S_III_APRP:       multi-state surface of C7H8S, 1-3 state,
                              availability: potential energy, gradient, nonadiabatic coupling vector 
                              functional form: APRP + PIP
                              corresponding surface in POTLIB: phsch3_III.f90
                              pecial emphsize on: C6H5SCH3 → C6H5S + CH3
                              ref: Y. Shu, and D. G. Truhlar, 
                                   "Improved potential energy surfaces of thioanisole and the 
                                   effect of upper surface variantions on the production distribution 
                                   upon photodissociation",
                                   Chem. Phys. 515, 737-743 (2018).
    ============================================================================================
    3.  C7H8S_S22_APRP:       multi-state surface of C7H8S, 1-3 state,
                              availability: potential energy, gradient, nonadiabatic coupling vector 
                              functional form: APRP + PIP
                              corresponding surface in POTLIB: phsch3_S-2.2.f90
                              pecial emphsize on: C6H5SCH3 → C6H5S + CH3
                              ref: Y. Shu, and D. G. Truhlar, 
                                   "Improved potential energy surfaces of thioanisole and the 
                                   effect of upper surface variantions on the production distribution 
                                   upon photodissociation",
                                   Chem. Phys. 515, 737-743 (2018).


    ==VERY IMPORTANT NOTICE==
    The atomic ordering specified in geom is required for C7H8S_APRP
         H10      H9   H14,15
          \      /      \ 
           C5---C4       C13---H16
          /      \      /
    H11--C6       C3---S12 
          \      /      
           C1---C2      
          /      \ 
         H7       H8

        """)

    return


def check(system, surface, geom):

    n_carbon=0
    n_hydrogen=0
    n_sulfer=0
    natoms=len(geom)

    if natoms!=16:
        print("number of atoms not equal 16 for C7H8S system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='S'):
                n_sulfer=n_sulfer+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
        if n_sulfer!=1:
            print("number of Sulfer atoms not equal 1 for C7H8S system")
            sys.exit()
        if n_hydrogen!=8:
            print("number of Hydrogen atoms not equal 8 for C7H8S system")
            sys.exit()
        if n_carbon!=7:
            print("number of Carbon atoms not equal 7 for C7H8S system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for C7H8S_APRP: input Cartesian should in order of CCCCCCHHHHHSCHHH
    if surface=='C7H8S_APRP' or surface=='C7H8S_III_APRP' or surface=='C7H8S_S22_APRP' or surface=='C7H8S_APRP_DPEM' or surface=='C7H8S_III_APRP_DPEM' or surface=='C7H8S_S22_APRP_DPEM':
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
