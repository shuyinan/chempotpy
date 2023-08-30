'''
F2H
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for F2H:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   F2H_LEPS_JOT_1972:       single-state surface of F2H, ground state,      P/G
    2.   F2H_LEPS_KNPRY_1966:     single-state surface of F2H, ground state,      P/G

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for F2H:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   F2H_LEPS_JOT_1972:       single-state surface of F2H, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: hf2jotii2001.f
                                   ref: N. Jonathon, S. Okuda, and D. Timlin,
                                        "Initial vibrational energy distributions determined 
                                        by infra-red chemiluminescence",
                                        Mol. Phys. 24, 1143-1164 (1972).
    ============================================================================================
    2.   F2H_LEPS_KNPRY_1966:      single-state surface of F2H, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: hcl2eleps2001.f
                                   ref: P. J. Kuntz, E. M. Nemth, J. C. Polanyi, 
                                        S. D. Rosner, and C. E. Young,
                                        "Energy Distribution Among Products of Exothermic 
                                        Reactions. II. Repulsive, Mixed, and Attractive 
                                        Energy Release", 
                                        J. Chem. Phys. 44, 1168-1184 (1966).

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_fluorine=0
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
        if n_hydrogen!=1:
            print("number of Hydrogen atoms not equal 1 for F2H system")
            sys.exit()
        if n_fluorine!=2:
            print("number of Fluorine atoms not equal 2 for F2H system")
            sys.exit()

    xyz=np.zeros((natoms,3))
 
    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='F2H_LEPS_JOT_1972' or surface=='F2H_LEPS_KNPRY_1966':
        geom_ordered=sorted(geom, key=lambda x: ("H", "F").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "F").index(x[0]))
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
