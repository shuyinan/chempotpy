'''
H3
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for H3:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces: 
    1.   H3_GEN_A2_2002:       single-state surface of H3, ground state,    P
    2.   H3_GEN_A3_2002:       single-state surface of H3, ground state,    P
    3.   H3_GEN_A4_2002:       single-state surface of H3, ground state,    P
    4.   H3_GEN_CCI_2002:      single-state surface of H3, ground state,    P
    5.   H3_GEN_BH_2009:       single-state surface of H3, ground state,    P
    6.   H3_GEN_BKMP3_1996:    single-state surface of H3, ground state,    P
    7.   H3_GEN_BKMP_1991:     single-state surface of H3, ground state,    P/G
    8.   H3_GEN_PK2_1964:      single-state surface of H3, ground state,    P/G
    9.   H3_GEN_TK_1984:       single-state surface of H3, ground state,    P/G

    Multi-State Surfaces:
    1.   H3_GEN_DMBE_1987:     multi-state surface of H3, 1-2 state,        P/G/D

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for O2H:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   H3_GEN_A2_2002:       single-state surface of H3, ground state,
                               availability: potential energy
                               functional form: GEN
                               corresponding surface in POTLIB: a2.f
                               ref: S. L. Mielke, B. C. Garrett, and K. A. Peterson,
                                    "A hierarchical family of global analytic Born–
                                    Oppenheimer potential energy surfaces for the H + H2 
                                    reaction ranging in quality from double-zeta to the 
                                    complete basis set limit",
                                    J. Chem. Phys. 116, 4142-4161 (2002).
    ============================================================================================
    2.   H3_GEN_A3_2002:       single-state surface of H3, ground state,
                               availability: potential energy
                               functional form: GEN
                               corresponding surface in POTLIB: a3.f
                               ref: S. L. Mielke, B. C. Garrett, and K. A. Peterson,
                                    "A hierarchical family of global analytic Born–
                                    Oppenheimer potential energy surfaces for the H + H2 
                                    reaction ranging in quality from double-zeta to the 
                                    complete basis set limit",
                                    J. Chem. Phys. 116, 4142-4161 (2002).
    ============================================================================================
    3.   H3_GEN_A4_2002:       single-state surface of H3, ground state,
                               availability: potential energy
                               functional form: GEN
                               corresponding surface in POTLIB: a4.f
                               ref: S. L. Mielke, B. C. Garrett, and K. A. Peterson,
                                    "A hierarchical family of global analytic Born–
                                    Oppenheimer potential energy surfaces for the H + H2 
                                    reaction ranging in quality from double-zeta to the 
                                    complete basis set limit",
                                    J. Chem. Phys. 116, 4142-4161 (2002).
    ============================================================================================
    4.   H3_GEN_CCI_2002:      single-state surface of H3, ground state,
                               availability: potential energy
                               functional form: GEN
                               corresponding surface in POTLIB: cci.f 
                               ref: S. L. Mielke, B. C. Garrett, and K. A. Peterson,
                                    "A hierarchical family of global analytic Born–
                                    Oppenheimer potential energy surfaces for the H + H2 
                                    reaction ranging in quality from double-zeta to the 
                                    complete basis set limit",
                                    J. Chem. Phys. 116, 4142-4161 (2002).
    ============================================================================================
    5.   H3_GEN_BH_2009:       single-state surface of H3, ground state,
                               availability: potential energy
                               functional form: GEN
                               corresponding surface in POTLIB: h3_bh.f 
                               ref: S. L. Mielke, D. W. Schwenke, G. C. Schatz, B. C. Garrett, 
                                    and K. A. Peterson,
                                    "Functional Representation for the Born−Oppenheimer Diagonal 
                                    Correction and Born−Huang Adiabatic Potential Energy Surfaces 
                                    for Isotopomers of H3",
                                    J. Phys. Chem. A 113, 4479-4488 (2009).
    ============================================================================================
    6.   H3_GEN_BKMP3_1996:    single-state surface of H3, ground state,
                               availability: potential energy
                               functional form: GEN
                               corresponding surface in POTLIB: h3bkmp3.f 
                               ref: A. I. Boothroyd, W. J. Keogh, P. G. Martin, and M. R. Peterson,
                                    "A refined H3 potential energy surface",
                                    J. Chem. Phys. 104, 7139-7152 (1996).
    ============================================================================================
    7.   H3_GEN_BKMP_1991:     single-state surface of H3, ground state, 
                               availability: potential energy
                               functional form: GEN
                               corresponding surface in POTLIB: h3bkmp.f
                               ref: A. I. Boothroyd, W. J. Keogh, P. G. Martin, and M. R. Peterson,
                                    "An improved H3 potential energy surface", 
                                    J. Chem. Phys. 95, 4343-4359 (1991).
    ========================================================================================
    8.   H3_GEN_PK2_1964:      single-state surface of H3, ground state,
                               availability: potential energy
                               functional form: GEN
                               corresponding surface in POTLIB: h3pk2-2001.f
                               ref: R. N. Porter, M. Karplus, 
                                    "Potential Energy Surface for H3",
                                    J. Chem. Phys. 40, 1105-1115 (1964).
    ========================================================================================
    9.   H3_GEN_TK_1984:       single-state surface of H3, ground state,
                               availability: potential energy
                               functional form: GEN
                               corresponding surface in POTLIB: h3tk2001.f
                               ref: D. G. Truhlar, and A. Kuppermann,
                                    "Exact and Approximate Quantum Mechanical Reaction Probabilities 
                                    and Rate Constants for the Collinear H + H2 Reaction",
                                    J. Chem. Phys. 56, 2232-2252 (1972).


    Multi-State Surfaces:
    1.   H3_GEN_DMBE_1987:     multi-state surface of H3, 1-2 state,        P/G/D
                               availability: potential energy
                               functional form: GEN
                               corresponding surface in POTLIB: h3dmbe2001.f 
                               ref: A. J. C. Varandas, F. B. Brown, C. A. Mead, D. G. Truhlar, 
                                    and N. C. Blais,
                                    "A double many‐body expansion of the two lowest‐energy potential 
                                    surfaces and nonadiabatic coupling for H3",
                                    J. Chem. Phys. 86, 6258-6269 (1987).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    natoms=len(geom)

    if natoms!=3:
        check=0
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
        if n_hydrogen!=3:
            print("number of Hydrogen atoms not equal 3 for H3 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    for iatom in range(natoms):
        for idir in range(3):
            xyz[iatom][idir]=geom[iatom][idir+1] 

    return xyz, rank

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
