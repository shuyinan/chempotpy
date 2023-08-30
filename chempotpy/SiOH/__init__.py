'''
SiOH
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for SiOH:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   SiOH_GEN_1988:     single-state surface of SiOH, ground state,    P/G

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for SiOH:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   SiOH_GEN_1988:     single-state surface of SiOH, ground state,
                            availability: potential energy, gradient
                            functional form: LEPS
                            corresponding surface in POTLIB: SiOH.f
                            ref: B. P. Feuston, S. H. Garofalini, 
                                 "Empirical three‐body potential for vitreous silica",
                                 J. Chem. Phys. 89, 5818-5824 (1988).
                                 D. A. Litton, S. H. Garofalini, 
                                 "Modeling of hydrophilic wafer bonding by molecular 
                                 dynamics simulations",
                                 J. Appl. Phys. 89, 6013-6023 (2001).

        """)

    return

def check(system, surface, geom):

    n_silicon=0
    n_hydrogen=0
    n_oxygen=0
    natoms=len(geom)

    if natoms==1:
        print("number of atoms not equal 1")
        sys.exit()
    elif natoms>10000:
        print("number of atoms not exceed 10000")
        sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='SiOH_GEN_1988':
        geom_ordered=sorted(geom, key=lambda x: ("Si", "O", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Si", "O", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    return xyz, rank_ordered



def nSinOnH(system, surface, geom):

    n_silicon=0
    n_hydrogen=0
    n_oxygen=0
    natoms=len(geom)

    if natoms==1:
        print("number of atoms not equal 1")
        sys.exit()
    elif natoms>10000:
        print("number of atoms not exceed 10000")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Si'):
                 n_silicon=n_silicon+1
            elif (geom[iatom][0]=='O'):
                 n_oxygen=n_oxygen+1

    return n_silicon, n_oxygen, n_hydrogen

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
