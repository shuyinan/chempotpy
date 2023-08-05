'''
MXH
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for MXH:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Multi-State Surfaces:
    1.   MXHWL_LEPS_ModelSurface:       multi-state surface of MXH, 1-2 states,    P/G/D
    2.   MXHSL_LEPS_ModelSurface:       multi-state surface of MXH, 1-2 states,    P/G/D
    3.   MXHSB_LEPS_ModelSurface:       multi-state surface of MXH, 1-2 states,    P/G/D
    4.   MXHWL_LEPS_ModelSurface_DPEM   multi-state surface of MXH, 1-2 states,    U/UG
    5.   MXHSL_LEPS_ModelSurface_DPEM   multi-state surface of MXH, 1-2 states,    U/UG
    6.   MXHSB_LEPS_ModelSurface_DPEM   multi-state surface of MXH, 1-2 states,    U/UG

    ==VERY IMPORTANT NOTICE==
    MXH is a model system, it is used to model a M atom interacts with diatomic XH molecule
    The MXH systems feature avoided crossings of the Landau-Zener type in which the two 
    diabatic potential energy surfaces cross with nonzero diabatic coupling. 

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
    1.   MXHWL_LEPS_ModelSurface:  multi-state surface of MXH, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: MXH_WLd.f
                                   ref: Y. L. Volobuev, M. D. Hack, M. S. Topaler, and D. G. Truhlar,
                                        "Continuous surface switching: An improved time-dependent 
                                        self-consistent-field method for nonadiabatic dynamics",
                                        J. Chem. Phys. 112, 9716-9726 (2000).
    ============================================================================================
    2.   MXHSL_LEPS_ModelSurface:  multi-state surface of MXH, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: MXH_SLd.f 
                                   ref: Y. L. Volobuev, M. D. Hack, M. S. Topaler, and D. G. Truhlar,
                                        "Continuous surface switching: An improved time-dependent 
                                        self-consistent-field method for nonadiabatic dynamics",
                                        J. Chem. Phys. 112, 9716-9726 (2000).
    ============================================================================================
    3.   MXHSB_LEPS_ModelSurface:  multi-state surface of MXH, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: MXH_SBd.f
                                   ref: Y. L. Volobuev, M. D. Hack, M. S. Topaler, and D. G. Truhlar,
                                        "Continuous surface switching: An improved time-dependent 
                                        self-consistent-field method for nonadiabatic dynamics",
                                        J. Chem. Phys. 112, 9716-9726 (2000).

    ==VERY IMPORTANT NOTICE==
    MXH is a model system, it is used to model a M atom interacts with diatomic XH molecule
    The MXH systems feature avoided crossings of the Landau-Zener type in which the two 
    diabatic potential energy surfaces cross with nonzero diabatic coupling. 

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_m=0
    n_x=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='M'):
                 n_m=n_m+1
            elif (geom[iatom][0]=='X'):
                 n_x=n_x+1

        if n_hydrogen!=1:
            print("number of Hydrogen atoms not equal 1 for MXH system")
            sys.exit()
        if n_m!=1:
            print("number of M atoms not equal 1 for MXH system")
            sys.exit()
        if n_x!=1:
            print("number of X atoms not equal 1 for MXH system")
            sys.exit()

    xyz=np.zeros((natoms,3))
   
    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='MXHWL_LEPS_ModelSurface' or surface=='MXHSL_LEPS_ModelSurface' or surface=='MXHSB_LEPS_ModelSurface':
        geom_ordered=sorted(geom, key=lambda x: ("M", "H", "X").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("M", "H", "X").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    if surface=='MXHWL_LEPS_ModelSurface_DPEM' or surface=='MXHSL_LEPS_ModelSurface_DPEM' or surface=='MXHSB_LEPS_ModelSurface_DPEM':
        geom_ordered=sorted(geom, key=lambda x: ("M", "H", "X").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("M", "H", "X").index(x[0]))
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
