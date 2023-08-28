'''
C2H6H
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for C2H6H:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  C2H6H_VBMM:           single-state surface of C2H6H, ground state,      P

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for C2H6H:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  C2H6H_VBMM:           single-state surface of C2H6H, ground state,
                              availability: potential energy
                              functional form: VBMM
                              corresponding surface in POTLIB: C2H6H-rev2019.f
                              special emphsize on: H + C2H6 reaction
                              ref: J. Espinosa-Garcia, M. Garcia-Chamorro, J. C. Corchado,
                                   S. Bhowmick, and Y. V. Suleimanov,
                                   "VTST and RPMD kinetics study of the nine-body X + C2H6 
                                   (X ≡ H, Cl, F) reactions based on analytical potential 
                                   energy surfaces"
                                   Phys. Chem. Chem. Phys. 22, 13790-13801 (2020).

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_carbon=0
    natoms=len(geom)

    if natoms!=9:
        print("number of atoms not equal 9 for C2H6H system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
        if n_carbon!=2:
            print("number of Carbon atoms not equal 2 for C2H6H system")
            sys.exit()
        if n_hydrogen!=7:
            print("number of Hydrogen atoms not equal 7 for C2H6H system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for C2H6H_VBMM: input Cartesian should in order of CCHHHHHHH
    if surface=='C2H6H_VBMM':
        geom_ordered=sorted(geom, key=lambda x: ("C", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("C", "H").index(x[0]))
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
