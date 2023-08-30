'''
AlmHn
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for AlmHn:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   AlmHn_VBO1:                single-state surface of AlmHn, ground state, P
    2.   AlmHn_VBO2:                single-state surface of AlmHn, ground state, P

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for MXH:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Multi-State Surfaces:
    1.   AlmHn_VBO1:                single-state surface of AlmHn, ground state,
                                    availability: potential energy, gradient
                                    functional form: GEN (Valence bond order)
                                    corresponding surface in POTLIB: vbo-beta1.f
                                    ref: M. Zhao, M. A. Iron, R. Staszewski, N. E. Schultz, 
                                         R. Valero, and D. G. Truhlar,
                                         "Valence–Bond Order (VBO): A New Approach to Modeling 
                                         Reactive Potential Energy Surfaces for Complex Systems, 
                                         Materials, and Nanoparticles",
                                         J. Chem. Theory Comput. 5, 594-604 (2009).
    ============================================================================================
    2.   AlmHn_VBO2:                single-state surface of AlmHn, ground state,
                                    availability: potential energy, gradient
                                    functional form: GEN (Valence bond order)
                                    corresponding surface in POTLIB: vbo2-beta1.f
                                    ref: M. Zhao, M. A. Iron, R. Staszewski, N. E. Schultz, 
                                         R. Valero, and D. G. Truhlar,
                                         "Valence–Bond Order (VBO): A New Approach to Modeling 
                                         Reactive Potential Energy Surfaces for Complex Systems, 
                                         Materials, and Nanoparticles",
                                         J. Chem. Theory Comput. 5, 594-604 (2009).


        """)

    return

def check(system, surface, geom):

    n_aluminium=0
    n_hydrogen=0
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

    geom_ordered=sorted(geom, key=lambda x: ("Al", "H").index(x[0]))
    rank_ordered=sorted(rank, key=lambda x: ("Al", "H").index(x[0]))
    for iatom in range(natoms):
        for idir in range(3):
            xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    return xyz, rank_ordered

def nAlnH(system, surface, geom):

    n_aluminium=0
    n_hydrogen=0
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
            elif (geom[iatom][0]=='Al'):
                 n_aluminium=n_aluminium+1

    return n_hydrogen, n_aluminium


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
