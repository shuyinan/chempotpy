'''
YRH
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for YRH:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Multi-State Surfaces:
    1.   YRH_LEPS_ModelSurface:        multi-state surface of YRH, 1-2 states,    P/G/D
    2.   YRH_LEPS_ModelSurface_DPEM:   multi-state surface of YRH, 1-2 states,    U/UG

    ==VERY IMPORTANT NOTICE==
    YRH is a model system, it is used to model a Y atom interacts with diatomic RH molecule
    The YRH systems feature interactions of the Rosen-Zener-Demkov type, by which we mean 
    cases where the two diabatic potential energy surfaces do not cross and are weakly coupled.

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
    1.   YRH_LEPS_ModelSurface:    multi-state surface of YRH, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: yrh_der.f 
                                   ref: A. W. Jasper, M. D. Hack, and D. G. Truhlar,
                                        "The treatment of classically forbidden electronic transitions in 
                                        semiclassical trajectory surface hopping calculations",
                                        J. Chem. Phys. 115, 1804-1816 (2001).

    ==VERY IMPORTANT NOTICE==
    YRH is a model system, it is used to model a Y atom interacts with diatomic RH molecule
    The YRH systems feature interactions of the Rosen-Zener-Demkov type, by which we mean 
    cases where the two diabatic potential energy surfaces do not cross and are weakly coupled.

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_y=0
    n_r=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Y'):
                 n_y=n_y+1
            elif (geom[iatom][0]=='R'):
                 n_r=n_r+1

        if n_hydrogen!=1:
            print("number of Hydrogen atoms not equal 1 for YRH system")
            sys.exit()
        if n_y!=1:
            print("number of Y atoms not equal 1 for YRH system")
            sys.exit()
        if n_r!=1:
            print("number of R atoms not equal 1 for YRH system")
            sys.exit()

    xyz=np.zeros((natoms,3))
 
    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)   

    if surface=='YRH_LEPS_ModelSurface' or surface=='YRH_LEPS_ModelSurface_DPEM': 
        geom_ordered=sorted(geom, key=lambda x: ("Y", "R", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Y", "R", "H").index(x[0]))
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
