'''
MCH
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for MCH:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Multi-State Surfaces:
    1.   MCHWL_LEPS_ModelSurface:      multi-state surface of MCH, 1-2 states,    P/G/D
    2.   MCHWB_LEPS_ModelSurface:      multi-state surface of MCH, 1-2 states,    P/G/D
    3.   MCHSL_LEPS_ModelSurface:      multi-state surface of MCH, 1-2 states,    P/G/D
    4.   MCHSB_LEPS_ModelSurface:      multi-state surface of MCH, 1-2 states,    P/G/D
    5.   MCHTL_LEPS_ModelSurface:      multi-state surface of MCH, 1-2 states,    P/G/D
    6.   MCHWL_LEPS_ModelSurface_DPEM: multi-state surface of MCH, 1-2 states,    U/UG
    7.   MCHWB_LEPS_ModelSurface_DPEM: multi-state surface of MCH, 1-2 states,    U/UG
    8.   MCHSL_LEPS_ModelSurface_DPEM: multi-state surface of MCH, 1-2 states,    U/UG
    9.   MCHSB_LEPS_ModelSurface_DPEM: multi-state surface of MCH, 1-2 states,    U/UG
    10.  MCHTL_LEPS_ModelSurface_DPEM: multi-state surface of MCH, 1-2 states,    U/UG


    ==VERY IMPORTANT NOTICE==
    MCH is a model system, it is used to model a M atom interacts with diatomic CH molecule
    The MCH systems feature conical intersections in which the two diabatic potential energy 
    surfaces cross at some geometries where the diabatic coupling is zero. 

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for MCH:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Multi-State Surfaces:
    1.   MCHWL_LEPS_ModelSurface:  multi-state surface of MCH, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: MCH_WL.f
                                   ref: A. W. Jasper, and D. G. Truhlar,
                                        "Conical intersections and semiclassical trajectories: Comparison to 
                                        accurate quantum dynamics and analyses of the trajectories",
                                        J. Chem. Phys. 122 044101 (2005).
    ============================================================================================
    2.   MCHWB_LEPS_ModelSurface:  multi-state surface of MCH, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: MCH_WB.f
                                   ref: A. W. Jasper, and D. G. Truhlar,
                                        "Conical intersections and semiclassical trajectories: Comparison to 
                                        accurate quantum dynamics and analyses of the trajectories",
                                        J. Chem. Phys. 122 044101 (2005).
    ============================================================================================
    3.   MCHSL_LEPS_ModelSurface:  multi-state surface of MCH, 1-2 states, 
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: MCH_SL.f
                                   ref: A. W. Jasper, and D. G. Truhlar,
                                        "Conical intersections and semiclassical trajectories: Comparison to 
                                        accurate quantum dynamics and analyses of the trajectories",
                                        J. Chem. Phys. 122 044101 (2005).
    ============================================================================================
    4.   MCHSB_LEPS_ModelSurface:  multi-state surface of MCH, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: MCH_SB.f
                                   ref: A. W. Jasper, and D. G. Truhlar,
                                        "Conical intersections and semiclassical trajectories: Comparison to 
                                        accurate quantum dynamics and analyses of the trajectories",
                                        J. Chem. Phys. 122 044101 (2005).
    ============================================================================================
    5.   MCHTL_LEPS_ModelSurface:  multi-state surface of MCH, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: MCH_TL.f
                                   ref: A. W. Jasper, and D. G. Truhlar,
                                        "Conical intersections and semiclassical trajectories: Comparison to 
                                        accurate quantum dynamics and analyses of the trajectories",
                                        J. Chem. Phys. 122 044101 (2005).

    ==VERY IMPORTANT NOTICE==
    MCH is a model system, it is used to model a M atom interacts with diatomic CH molecule
    The MCH systems feature conical intersections in which the two diabatic potential energy 
    surfaces cross at some geometries where the diabatic coupling is zero. 

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_m=0
    n_c=0
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
            elif (geom[iatom][0]=='C'):
                 n_c=n_c+1

        if n_hydrogen!=1:
            print("number of Hydrogen atoms not equal 1 for MCH system")
            sys.exit()
        if n_m!=1:
            print("number of M atoms not equal 1 for MCH system")
            sys.exit()
        if n_c!=1:
            print("number of C atoms not equal 1 for MCH system")
            sys.exit()

    xyz=np.zeros((natoms,3))
  
    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='MCHWL_LEPS_ModelSurface' or surface=='MCHWB_LEPS_ModelSurface' or surface=='MCHSL_LEPS_ModelSurface' or surface=='MCHSB_LEPS_ModelSurface' or surface=='MCHTL_LEPS_ModelSurface':
        geom_ordered=sorted(geom, key=lambda x: ("M", "H", "C").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("M", "H", "C").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    if surface=='MCHWL_LEPS_ModelSurface_DPEM' or surface=='MCHWB_LEPS_ModelSurface_DPEM' or surface=='MCHSL_LEPS_ModelSurface_DPEM' or surface=='MCHSB_LEPS_ModelSurface_DPEM' or surface=='MCHTL_LEPS_ModelSurface_DPEM':
        geom_ordered=sorted(geom, key=lambda x: ("M", "H", "C").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("M", "H", "C").index(x[0]))
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
