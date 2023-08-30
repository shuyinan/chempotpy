'''
OH2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for OH2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces: 
    1.   OH2_Ap_GEN_M21986:    single-state surface of OH2, A' symmetry, ground state,    P/G
    2.   OH2_App_GEN_M21986:   single-state surface of OH2, A'' symmetry, ground state,   P/G
    3.   OH2_App_GEN_M3n1988:  single-state surface of OH2, A'' symmetry, ground state,   P/G
    4.   OH2_App_GEN_M3na1988: single-state surface of OH2, A'' symmetry, ground state,   P/G
    5.   OH2_GEN_DIM1981:      single-state surface of OH2, ground state,                 P/G
    6.   OH2_LEPS_JW1977:      single-state surface of OH2, ground state,                 P/G
    7.   OH2_LEPS_JWS1985:     single-state surface of OH2, ground state,                 P/G
    8.   OH2_GEN_POL1980:      single-state surface of OH2, ground state,                 P/G
    9.   OH2_GEN_SL1979:       single-state surface of OH2, ground state,                 P/G
 
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
    1.   OH2_Ap_GEN_M21986:    single-state surface of OH2, A' symmetry, ground state,
                               availability: potential energy, gradient
                               functional form: GEN
                               corresponding surface in POTLIB: oh2m2ap2001.f
                               ref: B. C. Garrett, and D. G. Truhlar,
                                    "Thermal and state‐selected rate constant calculations 
                                    for O(3p) + H2 → OH + H and isotopic analogs",
                                    Int. J. Quantum Chem. 29, 1463-1482 (1986).  
    ============================================================================================
    2.   OH2_App_GEN_M21986:   single-state surface of OH2, A'' symmetry, ground state,
                               availability: potential energy, gradient
                               functional form: GEN
                               corresponding surface in POTLIB: oh2m2adp2001.f
                               ref: B. C. Garrett, and D. G. Truhlar,
                                    "Thermal and state‐selected rate constant calculations 
                                    for O(3p) + H2 → OH + H and isotopic analogs",
                                    Int. J. Quantum Chem. 29, 1463-1482 (1986). 
    ============================================================================================
    3.   OH2_App_GEN_M3n1988:  single-state surface of OH2, A'' symmetry, ground state,
                               availability: potential energy, gradient
                               functional form: GEN
                               corresponding surface in POTLIB: oh2m3n2001.f
                               ref: T. Joseph, D. G. Truhlar, and B. C. Garrett,
                                    "Improved potential energy surfaces for the reaction O(3P)
                                    + H2 → OH + H", 
                                    J. Chem. Phys. 88, 6982-6990 (1988).
    ============================================================================================
    4.   OH2_App_GEN_M3na1988: single-state surface of OH2, A'' symmetry, ground state,
                               availability: potential energy, gradient
                               functional form: GEN
                               corresponding surface in POTLIB: oh2m3n2001.f
                               ref: T. Joseph, D. G. Truhlar, and B. C. Garrett,
                                    "Improved potential energy surfaces for the reaction O(3P)
                                    + H2 → OH + H", 
                                    J. Chem. Phys. 88, 6982-6990 (1988).
    ============================================================================================
    5.   OH2_GEN_DIM1981:      single-state surface of OH2, ground state,
                               availability: potential energy, gradient
                               functional form: GEN
                               corresponding surface in POTLIB: oh2dim2001.f
                               ref: G. C. Schatz, A. F. Wagner, S. P. Walch, J. M. Bowman, 
                                    "A comparative study of the reaction dynamics of several 
                                    potential energy surfaces of O(3P) + H2 → OH + H. I",
                                    J. Chem. Phys. 74, 4984-4996 (1981).
    ============================================================================================
    6.   OH2_LEPS_JW1977:      single-state surface of OH2, ground state,
                               availability: potential energy, gradient
                               functional form: GEN
                               corresponding surface in POTLIB: oh2jw2001.f
                               ref: B. R. Johnson, and N. W. Winter,
                                    "Classical trajectory study of the effect of vibrational 
                                    energy on the reaction of molecular hydrogen with atomic 
                                    oxygen",
                                    J. Chem. Phys. 66, 4116-4120 (1977).
    ============================================================================================
    7.   OH2_LEPS_JWS1985:     single-state surface of OH2, ground state,
                               availability: potential energy, gradient
                               functional form: GEN
                               corresponding surface in POTLIB: oh2jws2001.f
                               ref: G. C. Schatz,
                                    "A coupled states distorted wave study of the O(3P) + H2 
                                    (D2, HD, DH) reaction", 
                                    J. Chem. Phys. 83, 5677-5686 (1985).
    ============================================================================================
    8.   OH2_GEN_POL1980:      single-state surface of OH2, ground state,
                               availability: potential energy, gradient
                               functional form: GEN
                               corresponding surface in POTLIB: oh2pol2001.f
                               ref: S. P. Walch, T. H. Dunning, Jr., F. W. Bobrowicz, R. Raffenetti,
                                    "A theoretical study of the potential energy surface for O(3P) + H2",
                                    J. Chem. Phys. 72, 406-415 (1980).
                                    A. F. Wagner, G. C. Schatz, and T. H. Dunning, Jr.,
                                    "Theoretical studies of the O + H2 reaction",
                                    J. Chem. Phys. 72, 2894-2896 (1980).
                                    A. F. Wagner, G. C. Schatz, J. M. Bowman, 
                                    "The evaluation of fitting functions for the representation of an 
                                    O(3P) + H2 potential energy surface. I ",
                                    J. Chem. Phys. 74, 4960-4983 (1981).
    ============================================================================================
    9.   OH2_GEN_SL1979:       single-state surface of OH2, ground state, 
                               availability: potential energy, gradient
                               functional form: GEN
                               corresponding surface in POTLIB: oh2sl2001.f
                               ref: R. Schinke and W. A. Lester, Jr.,
                                    "Trajectory study of O + H2 reactions on fitted ab initio surfaces I: 
                                    Triplet case",
                                    J. Chem. Phys. 70, 4893-4902 (1979).


        """)

    return

def check(system, surface, geom):

    n_oxygen=0
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
            elif (geom[iatom][0]=='O'):
                 n_oxygen=n_oxygen+1
        if n_hydrogen!=2:
            print("number of Hydrogen atoms not equal 2 for OH2 system")
            sys.exit()
        if n_oxygen!=1:
            print("number of Oxygen atoms not equal 1 for OH2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='OH2_Ap_GEN_M21986' or surface=='OH2_App_GEN_M21986' or surface=='OH2_App_GEN_M3n1988' or surface=='OH2_App_GEN_M3na1988' or surface=='OH2_GEN_DIM1981' or surface=='OH2_LEPS_JW1977' or surface=='OH2_LEPS_JWS1985' or surface=='OH2_GEN_POL1980' or surface=='OH2_GEN_SL1979':
        geom_ordered=sorted(geom, key=lambda x: ("O", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("O", "H").index(x[0]))
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
