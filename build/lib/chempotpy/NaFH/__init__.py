'''
NaFH
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for NaFH:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Multi-State Surfaces:
    1.   NaFH_LEPS_JHCTP_2001:       multi-state surface of NaFH, 1-2 states,    P/G/D
    2.   NaFH_LEPS_TTYPP1_1998:      multi-state surface of NaFH, 1-2 states,    P
    3.   NaFH_LEPS_TTYPP2_1998:      multi-state surface of NaFH, 1-2 states,    P/G/D 
    4.   NaFH_LEPS_JHCTP_2001_DPEM:  multi-state surface of NaFH, 1-2 states,    U/UG
    5.   NaFH_LEPS_TTYPP1_1998_DPEM: multi-state surface of NaFH, 1-2 states,    U/UG
    6.   NaFH_LEPS_TTYPP2_1998_DPEM: multi-state surface of NaFH, 1-2 states,    U/UG
 
        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for BrH2:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Multi-State Surfaces:
    1.   NaFH_LEPS_JHCTP_2001:     multi-state surface of NaFH, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: nafh2d.f 
                                   ref: A. W. Jasper, M. D. Hack, A. Chakraborty, D. G. Truhlar, and P. Piecuch,
                                        "Photodissociation of LiFH and NaFH van der Waals complexes: 
                                        A semiclassical trajectory study",
                                        J. Chem. Phys. 115, 7945-7952 (2001).
    ============================================================================================
    2.   NaFH_LEPS_TTYPP1_1998:    multi-state surface of NaFH, 1-2 states,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: nafh1v.f
                                   ref: M. S. Topaler, D. G. Truhlar, X. Y. Chang, P. Piecuch, and J. C. Polanyi,
                                        "Potential energy surfaces of NaFH",
                                        J. Chem. Phys. 108, 5349-5377 (1998).
    ============================================================================================
    3.   NaFH_LEPS_TTYPP2_1998:    multi-state surface of NaFH, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: nafh1v.f
                                   ref: M. S. Topaler, D. G. Truhlar, X. Y. Chang, P. Piecuch, and J. C. Polanyi,
                                        "Potential energy surfaces of NaFH",
                                        J. Chem. Phys. 108, 5349-5377 (1998).

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_fluorine=0
    n_sodium=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='F'):
                 n_fluorine=n_fluorine+1
            elif (geom[iatom][0]=='Na'):
                 n_sodium=n_sodium+1

        if n_hydrogen!=1:
            print("number of Hydrogen atoms not equal 1 for NaFH system")
            sys.exit()
        if n_fluorine!=1:
            print("number of Fluorine atoms not equal 1 for NaFH system")
            sys.exit()
        if n_sodium!=1:
            print("number of Sodium atoms not equal 1 for NaFH system")
            sys.exit()

    xyz=np.zeros((natoms,3))
     
    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='NaFH_LEPS_JHCTP_2001' or surface=='NaFH_LEPS_TTYPP1_1998' or surface=='NaFH_LEPS_TTYPP2_1998' or surface=='NaFH_LEPS_JHCTP_2001_DPEM' or surface=='NaFH_LEPS_TTYPP1_1998_DPEM' or surface=='NaFH_LEPS_TTYPP2_1998_DPEM':
        geom_ordered=sorted(geom, key=lambda x: ("Na", "H", "F").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Na", "H", "F").index(x[0]))
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
