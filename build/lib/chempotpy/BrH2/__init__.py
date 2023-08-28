'''
BrH2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for BrH2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   BrH2_LEPS_LTBZ1_1995:       single-state surface of BrH2, ground state,      P/G
    2.   BrH2_LEPS_W_1959:           single-state surface of BrH2, ground state,      P/G
    3.   BrH2_LEPS_C_1982:           single-state surface of BrH2, ground state,      P/G

    Multi-State Surfaces:
    1.   BrH2_LEPS_LTBZSO_1995:      multi-state surface of BrH2, 1-2 states,         P
    2.   BrH2_LEPS_LTBZ2_1995:       multi-state surface of BrH2, 1-2 states,         P/G/D
    3.   BrH2_LEPS_LTBZSO_1995_DPEM: multi-state surface of BrH2, 1-2 states,         U
    4.   BrH2_LEPS_LTBZ2_1995_DPEM:  multi-state surface of BrH2, 1-2 states,         U/UG 

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

    Single-State Surfaces:
    1.   BrH2_LEPS_LTBZ1_1995:     single-state surface of BrH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: brh2sec2001.f 
                                   ref: G. C. Lynch, D. G. Truhlar, F. B. Brown, and
                                        J.-g. Zhao, 
                                        "A New Potential Energy Surface for H2Br and Its 
                                        Use To Calculate Branching Ratios and Kinetic Isotope 
                                        Effects for the H + HBr Reaction",
                                        J. Phys. Chem. 99, 207-225 (1995).
    ============================================================================================
    2.   BrH2_LEPS_W_1959:         single-state surface of BrH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: brh2leps.f
                                   ref: R. Weston,
                                        "H3 Activated Complex and the Rate of Reaction of Hydrogen 
                                        Atoms with Hydrogen Molecules",
                                        J. Chem. Phys. 31, 892-898 (1959).
    ============================================================================================
    3.   BrH2_LEPS_C_1982:         single-state surface of BrH2, ground state, 
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: brh2dim3c2001.f 
                                   ref: D. C. Clary,
                                        "Exchange reactions of hydrogen halides with hydrogenic 
                                        atoms",
                                        Chem. Phys. 71, 117-125 (1982).

    Multi-State Surfaces:
    1.   BrH2_LEPS_LTBZSO_1995:    single-state surface of BrH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: brh2so.f 
                                   ref: G. C. Lynch, D. G. Truhlar, F. B. Brown, and
                                        J.-g. Zhao, 
                                        "A New Potential Energy Surface for H2Br and Its 
                                        Use To Calculate Branching Ratios and Kinetic Isotope 
                                        Effects for the H + HBr Reaction",
                                        J. Phys. Chem. 99, 207-225 (1995).
    ============================================================================================
    2.   BrH2_LEPS_LTBZ2_1995:     multi-state surface of BrH2, 1-2 states,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: brh2secd.f 
                                   ref: G. C. Lynch, D. G. Truhlar, F. B. Brown, and
                                        J.-g. Zhao, 
                                        "A New Potential Energy Surface for H2Br and Its 
                                        Use To Calculate Branching Ratios and Kinetic Isotope 
                                        Effects for the H + HBr Reaction",
                                        J. Phys. Chem. 99, 207-225 (1995).

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_bromine=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Br'):
                 n_bromine=n_bromine+1
        if n_hydrogen!=2:
            print("number of Hydrogen atoms not equal 2 for BrH2 system")
            sys.exit()
        if n_bromine!=1:
            print("number of Bromine atoms not equal 1 for BrH2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))
     
    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='BrH2_LEPS_LTBZ1_1995' or surface=='BrH2_LEPS_W_1959' or surface=='BrH2_LEPS_C_1982':
        geom_ordered=sorted(geom, key=lambda x: ("Br", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Br", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]


    if surface=='BrH2_LEPS_LTBZSO_1995' or surface=='BrH2_LEPS_LTBZ2_1995' or surface=='BrH2_LEPS_LTBZSO_1995_DPEM' or surface=='BrH2_LEPS_LTBZ2_1995_DPEM':
        geom_ordered=sorted(geom, key=lambda x: ("Br", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Br", "H").index(x[0]))
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
