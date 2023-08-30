'''
ALC
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for ALC:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   ALC_LEPS_ModelSurface:     single-state surface of ALC, ground state, P/G
    2.   ALC5m0_LEPS_ModelSurface:  single-state surface of ALC, ground state, P/G
    3.   ALC6m0_LEPS_ModelSurface:  single-state surface of ALC, ground state, P/G
    4.   ALC7m0_LEPS_ModelSurface:  single-state surface of ALC, ground state, P/G
    5.   ALC7m1_LEPS_ModelSurface:  single-state surface of ALC, ground state, P/G
    6.   ALC7m4_LEPS_ModelSurface:  single-state surface of ALC, ground state, P/G
    7.   ALC7m5_LEPS_ModelSurface:  single-state surface of ALC, ground state, P/G
    8.   ALC7m7_LEPS_ModelSurface:  single-state surface of ALC, ground state, P/G
    9.   ALC7m10_LEPS_ModelSurface: single-state surface of ALC, ground state, P/G
    10.  ALC7p2_LEPS_ModelSurface:  single-state surface of ALC, ground state, P/G
    11.  ALC7p3_LEPS_ModelSurface:  single-state surface of ALC, ground state, P/G
    12.  ALC7p5_LEPS_ModelSurface:  single-state surface of ALC, ground state, P/G


    ==VERY IMPORTANT NOTICE==
    ALC is a model system for three-body model reactionL: A + L-C -> A-L + C, where A is 
    the acceptor group, L is the H or D, and C denotes the methyl donor. 

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for ALC:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Multi-State Surfaces:
    1.   ALC_LEPS_ModelSurface:     single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc.f
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).
    ============================================================================================
    2.   ALC5m0_LEPS_ModelSurface:  single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc5m0.f 
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).
    ============================================================================================
    3.   ALC6m0_LEPS_ModelSurface:  single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc6m0.f
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).
    ===========================================================================================
    4.   ALC7m0_LEPS_ModelSurface:  single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc7m0.f
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).
    ===========================================================================================
    5.   ALC7m1_LEPS_ModelSurface:  single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc7m1.f
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).
    ===========================================================================================
    6.   ALC7m4_LEPS_ModelSurface:  single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc7m4.f
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).
    ===========================================================================================
    7.   ALC7m5_LEPS_ModelSurface:  single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc7m5.f
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).
    ===========================================================================================
    8.   ALC7m7_LEPS_ModelSurface:  single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc7m7.f
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).
    ===========================================================================================
    9.   ALC7m10_LEPS_ModelSurface: single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc7m10.f 
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).
    ===========================================================================================
    10.  ALC7p2_LEPS_ModelSurface:  single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc7p2.f
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).
    ===========================================================================================
    11.  ALC7p3_LEPS_ModelSurface:  single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc7p3.f
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).
    ===========================================================================================
    12.  ALC7p5_LEPS_ModelSurface:  single-state surface of ALC, ground state,
                                    availability: potential energy, gradient
                                    functional form: LEPS
                                    corresponding surface in POTLIB: chc7p5.f
                                    ref: M. M. Kreevoy, D. Ostovic, D. G. Truhlar, and B. C. Garrett,
                                         "Phenomenological Manifestations of Large-Curvature Tunneling 
                                         In Hydride-Transfer Reactions",
                                         J. Phys. Chem. 90, 3766-3774 (1986).


    ==VERY IMPORTANT NOTICE==
    ALC is a model system for three-body model reaction: A + L-C -> A-L + C, where A is 
    the acceptor group, L is the H or D, and C denotes the methyl donor.


        """)

    return

def check(system, surface, geom):

    n_A=0
    n_L=0
    n_C=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='A'):
                 n_A=n_A+1
            elif (geom[iatom][0]=='L'):
                 n_L=n_L+1
            elif (geom[iatom][0]=='C'):
                 n_C=n_C+1

        if n_A!=1:
            print("number of acceptor group (A) is not equal 1 for ALC system")
            sys.exit()
        if n_L!=1:
            print("number of hydrogen/deuterium atom (L) is not equal 1 for ALC system")
            sys.exit()
        if n_C!=1:
            print("number of methyl donor (C) is not equal 1 for ALC system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='ALC_LEPS_ModelSurface' or surface=='ALC5m0_LEPS_ModelSurface' or surface=='ALC6m0_LEPS_ModelSurface' or surface=='ALC7m0_LEPS_ModelSurface' or surface=='ALC7m1_LEPS_ModelSurface' or surface=='ALC7m4_LEPS_ModelSurface' or surface=='ALC7m5_LEPS_ModelSurface' or surface=='ALC7m7_LEPS_ModelSurface' or surface=='ALC7m10_LEPS_ModelSurface' or surface=='ALC7p2_LEPS_ModelSurface' or surface=='ALC7p3_LEPS_ModelSurface' or surface=='ALC7p5_LEPS_ModelSurface': 
        geom_ordered=sorted(geom, key=lambda x: ("A", "L", "C").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("A", "L", "C").index(x[0]))
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
