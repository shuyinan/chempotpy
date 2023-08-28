'''
NH3
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for NH3:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Multi-State Surfaces:
    1.  NH3_GEN_FFW1_2006:      multi-state surface of NH3, 1-2 state,     P/G/D
    2.  NH3_GEN_FFW2_2007:      multi-state surface of NH3, 1-2 state,     P/G/D
    3.  NH3_GEN_FFW1_2006_DPEM: multi-state surface of NH3, 1-2 state,     U/UG
    4.  NH3_GEN_FFW2_2007_DPEM: multi-state surface of NH3, 1-2 state,     U/UG

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for NH3:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  NH3_GEN_FFW1_2006:    multi-state surface of NH3, 1-2 state,
                              availability: potential energy, gradient, nonadiabatic coupling vector
                              functional form: GEN
                              corresponding surface in POTLIB: nh3code1.f 
                              ref: S. Nangia, and D. G. Truhlar,
                                   "Direct calculation of coupled diabatic potential
                                   -energy surfaces for ammonia and mapping of a four
                                   -dimensional conical intersection seam",
                                   J. Chem. Phys. 124, 124309 (2006).
    ============================================================================================
    2.  NH3_GEN_FFW2_2007:    multi-state surface of NH3, 1-2 state,
                              availability: potential energy, gradient, nonadiabatic coupling vector
                              functional form: GEN
                              corresponding surface in POTLIB: nh3code2.f
                              ref: Z. H. Li, R. Valero, and D. G. Truhlar,
                                   "Improved direct diabatization and coupled potential 
                                   energy surfaces for the photodissociation of ammonia",
                                   Theor. Chem. Acc. 118, 9-24 (2007).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_nitrogen=0
    natoms=len(geom)

    if natoms!=4:
        print("number of atoms not equal 4 for NH3 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='N'):
                 n_nitrogen=n_nitrogen+1
        if n_hydrogen!=3:
            print("number of Hydrogen atoms not equal 3 for NH3 system")
            sys.exit()
        if n_nitrogen!=1:
            print("number of Nitrogen atoms not equal 1 for NH4 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for NH3_GEN_FFW1_2006 and NH3_GEN_FFW2_2006: input Cartesian should in order of NHHH
    if surface=='NH3_GEN_FFW1_2006' or surface=='NH3_GEN_FFW2_2007' or surface=='NH3_GEN_FFW1_2006_DPEM' or surface=='NH3_GEN_FFW2_2007_DPEM':
        geom_ordered=sorted(geom, key=lambda x: ("N", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("N", "H").index(x[0]))
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
