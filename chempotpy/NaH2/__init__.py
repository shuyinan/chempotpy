'''
NaH2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for NaH2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Multi-State Surfaces:
    1.   NaH2_LEPS_7S_2000:        multi-state surface of NaH2, 1-2 states,    P
    2.   NaH2_LEPS_7SD_2000:       multi-state surface of NaH2, 1-2 states,    P/G/D
    3.   NaH2_LEPS_7L_2000:        multi-state surface of NaH2, 1-2 states,    P
    4.   NaH2_LEPS_7LD_2000:       multi-state surface of NaH2, 1-2 states,    P/G/D
    5.   NaH2_LEPS_6_2000:         multi-state surface of NaH2, 1-2 states,    P/G/D
    6.   NaH2_LEPS_5G_1992:        multi-state surface of NaH2, 1-2 states,    P
    7.   NaH2_LEPS_5F_1992:        multi-state surface of NaH2, 1-2 states,    P
    8.   NaH2_LEPS_7S_2000_DPEM:   multi-state surface of NaH2, 1-2 states,    U
    9.   NaH2_LEPS_7SD_2000_DPEM:  multi-state surface of NaH2, 1-2 states,    U/UG
    10.  NaH2_LEPS_7L_2000_DPEM:   multi-state surface of NaH2, 1-2 states,    U
    11.  NaH2_LEPS_7LD_2000_DPEM:  multi-state surface of NaH2, 1-2 states,    U/UG
    12.  NaH2_LEPS_6_2000_DPEM:    multi-state surface of NaH2, 1-2 states,    U/UG
    13.  NaH2_LEPS_5G_1992_DPEM:   multi-state surface of NaH2, 1-2 states,    U
    14.  NaH2_LEPS_5F_1992_DPEM:   multi-state surface of NaH2, 1-2 states,    U


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
    1.   NaH2_LEPS_7S_2000:        multi-state surface of NaH2, 1-2 states,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: nah2_7s.f
                                   ref: M. D. Hack, A. W. Jasper, Y. L. Volobuev, D. W. Schwenke, 
                                        and D. G. Truhlar,
                                        "Do Semiclassical Trajectory Theories Provide an Accurate Picture of 
                                        Radiationless Decay for Systems with Accessible Surface Crossings?",
                                        J. Phy.s Chem. A 104, 217-232 (2000)
    ============================================================================================
    2.   NaH2_LEPS_7SD_2000:       multi-state surface of NaH2, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: nah2_7sd.f
                                   ref: M. D. Hack, A. W. Jasper, Y. L. Volobuev, D. W. Schwenke, 
                                        and D. G. Truhlar,
                                        "Do Semiclassical Trajectory Theories Provide an Accurate Picture of 
                                        Radiationless Decay for Systems with Accessible Surface Crossings?",
                                        J. Phy.s Chem. A 104, 217-232 (2000)
    ============================================================================================
    3.   NaH2_LEPS_7L_2000:        multi-state surface of NaH2, 1-2 states,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: nah2_7l.f
                                   ref: M. D. Hack, A. W. Jasper, Y. L. Volobuev, D. W. Schwenke, 
                                        and D. G. Truhlar,
                                        "Do Semiclassical Trajectory Theories Provide an Accurate Picture of 
                                        Radiationless Decay for Systems with Accessible Surface Crossings?",
                                        J. Phy.s Chem. A 104, 217-232 (2000)
    ============================================================================================
    4.   NaH2_LEPS_7LD_2000:       multi-state surface of NaH2, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: nah2_7ld.f
                                   ref: M. D. Hack, A. W. Jasper, Y. L. Volobuev, D. W. Schwenke, 
                                        and D. G. Truhlar,
                                        "Do Semiclassical Trajectory Theories Provide an Accurate Picture of 
                                        Radiationless Decay for Systems with Accessible Surface Crossings?",
                                        J. Phy.s Chem. A 104, 217-232 (2000)
    ============================================================================================
    5.   NaH2_LEPS_6_2000:         multi-state surface of NaH2, 1-2 states,
                                   availability: potential energy, gradient, and nonadiabatic coupling vector 
                                   functional form: LEPS
                                   corresponding surface in POTLIB: nah26_der.f
                                   ref: Y. L. Volobuev, M. D. Hack, M. S. Topaler, and D. G. Truhlar,
                                        "Continuous surface switching: An improved time-dependent 
                                        self-consistent-field method for nonadiabatic dynamics",
                                        J. Chem. Phys. 112, 9716-9726 (2000).
    ============================================================================================
    6.   NaH2_LEPS_5G_1992:        multi-state surface of NaH2, 1-2 states,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: nah25g.f 
                                   ref: P. Halvick, and D. G. Truhlar, 
                                        "A new diabatic representation of the coupled potential energy surfaces 
                                        for Na(3p 2P) + H2 → Na(3s 2S) + H2 or NaH + H",
                                        J. Chem. Phys. 96, 2895-2909 (1992); E 100, 4718 (1994).
    ============================================================================================
    7.   NaH2_LEPS_5F_1992:        multi-state surface of NaH2, 1-2 states,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: nah25f.f
                                   ref: P. Halvick, and D. G. Truhlar, 
                                        "A new diabatic representation of the coupled potential energy surfaces 
                                        for Na(3p 2P) + H2 → Na(3s 2S) + H2 or NaH + H",
                                        J. Chem. Phys. 96, 2895-2909 (1992); E 100, 4718 (1994).

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_sodium=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Na'):
                 n_sodium=n_sodium+1

        if n_hydrogen!=2:
            print("number of Hydrogen atoms not equal 2 for NaH2 system")
            sys.exit()
        if n_sodium!=1:
            print("number of Sodium atoms not equal 1 for NaH2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))
     
    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='NaH2_LEPS_7S_2000' or surface=='NaH2_LEPS_7SD_2000' or surface=='NaH2_LEPS_7L_2000' or surface=='NaH2_LEPS_7LD_2000' or surface=='NaH2_LEPS_6_2000' or surface=='NaH2_LEPS_5G_1992' or surface=='NaH2_LEPS_5F_1992':
        geom_ordered=sorted(geom, key=lambda x: ("Na", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Na", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    if surface=='NaH2_LEPS_7S_2000_DPEM' or surface=='NaH2_LEPS_7SD_2000_DPEM' or surface=='NaH2_LEPS_7L_2000_DPEM' or surface=='NaH2_LEPS_7LD_2000_DPEM' or surface=='NaH2_LEPS_6_2000_DPEM' or surface=='NaH2_LEPS_5G_1992_DPEM' or surface=='NaH2_LEPS_5F_1992_DPEM':
        geom_ordered=sorted(geom, key=lambda x: ("Na", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Na", "H").index(x[0]))
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
