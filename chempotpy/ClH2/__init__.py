'''
ClH2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for ClH2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   ClH2_LEPS_GTM1_1981:      single-state surface of ClH2, ground state,      P/G
    2.   ClH2_LEPS_GTM2_1981:      single-state surface of ClH2, ground state,      P/G
    3.   ClH2_LEPS_GTM3_1981:      single-state surface of ClH2, ground state,      P/G
    4.   ClH2_LEPS_GTM4_1982:      single-state surface of ClH2, ground state,      P/G
    5.   ClH2_LEPS_PK1_1966:       single-state surface of ClH2, ground state,      P/G
    6.   ClH2_LEPS_PK2_1966:       single-state surface of ClH2, ground state,      P/G
    7.   ClH2_LEPS_STSBLTG1_1989:  single-state surface of ClH2, ground state,      P/G
    8.   ClH2_LEPS_STSBLTG2_1989:  single-state surface of ClH2, ground state,      P/G
    9.   ClH2_LEPS_KNPRY_1966:     single-state surface of ClH2, ground state,      P/G
    10.  ClH2_LEPS_ALTG1_1996:     single-state surface of ClH2, ground state,      P/G
    11.  ClH2_LEPS_ALTG2_1996:     single-state surface of ClH2, ground state,      P
    12.  ClH2_LEPS_SPK_1973:       single-state surface of ClH2, ground state,      P/G
    13.  ClH2_LEPS_LB_1980:        single-state surface of ClH2, ground state,      P/G
    14.  ClH2_LEPS_TTGI1_1985:     single-state surface of ClH2, ground state,      P/G
    15.  ClH2_LEPS_TTGI2_1985:     single-state surface of ClH2, ground state,      P/G

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for ClH2:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   ClH2_LEPS_GTM1_1981:      single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2s1.f
                                   ref: B. C. Garrett, D. G. Truhlar, and A. W. Magnuson,
                                        "Variational transition state theory and vibrationally 
                                        adiabatic transmission coefficients for kinetic isotope
                                        effects in the Cl–H–H reaction system",
                                        J. Chem. Phys. 74, 1029-1043 (1981).
    ============================================================================================
    2.   ClH2_LEPS_GTM2_1981:      single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2s2.f
                                   ref: B. C. Garrett, D. G. Truhlar, and A. W. Magnuson,
                                        "Variational transition state theory and vibrationally 
                                        adiabatic transmission coefficients for kinetic isotope
                                        effects in the Cl–H–H reaction system",
                                        J. Chem. Phys. 74, 1029-1043 (1981).
    ============================================================================================
    3.   ClH2_LEPS_GTM3_1981:      single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2s2.f
                                   ref: B. C. Garrett, D. G. Truhlar, and A. W. Magnuson,
                                        "Variational transition state theory and vibrationally 
                                        adiabatic transmission coefficients for kinetic isotope
                                        effects in the Cl–H–H reaction system",
                                        J. Chem. Phys. 74, 1029-1043 (1981).
    ============================================================================================
    4.   ClH2_LEPS_GTM4_1982:      single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2alab.f
                                   ref: B. C. Garrett, D. G. Truhlar, and A. W. Magnuson,
                                        "New semiempirical method of modeling potential energy 
                                        surfaces for generalized TST and application to the 
                                        kinetic isotope effects in the Cl–H–H system",
                                        J. Chem. Phys. 76, 2321-2331 (1982).
    ============================================================================================
    5.   ClH2_LEPS_PK1_1966:       single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2si.f
                                   ref: A. Persky, and F. S. Klein,
                                        "Kinetic Isotope Effects in the Reaction between Atomic 
                                        Chlorine and Molecular Hydrogen. Tunnel Coefficients of the 
                                        Hydrogen Atom through an Asymmetric Potential Barrier",
                                        J. Chem. Phys. 44, 3617-3626 (1966).
    ============================================================================================
    6.   ClH2_LEPS_PK2_1966:       single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2sii.f
                                   ref: A. Persky, and F. S. Klein,
                                        "Kinetic Isotope Effects in the Reaction between Atomic 
                                        Chlorine and Molecular Hydrogen. Tunnel Coefficients of the 
                                        Hydrogen Atom through an Asymmetric Potential Barrier",
                                        J. Chem. Phys. 44, 3617-3626 (1966).
    ============================================================================================
    7.   ClH2_LEPS_STSBLTG1_1989:  single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2gq.f
                                   ref: D. W. Schwenke, S. C. Tucker, R. Steckler, F. B. Brown,
                                        G. C. Lynch, D. G. Truhlar, and B. C. Garrett,
                                        "Global potential‐energy surfaces for H2Cl", 
                                        J. Chem. Phys. 90, 3110-3120 (1989).
    ============================================================================================
    8.   ClH2_LEPS_STSBLTG2_1989:  single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2gqq.f
                                   ref: D. W. Schwenke, S. C. Tucker, R. Steckler, F. B. Brown,
                                        G. C. Lynch, D. G. Truhlar, and B. C. Garrett,
                                        "Global potential‐energy surfaces for H2Cl", 
                                        J. Chem. Phys. 90, 3110-3120 (1989).
    ============================================================================================
    9.   ClH2_LEPS_KNPRY_1966:     single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2g1.f 
                                   ref: P. J. Kuntz, E. M. Nemth, J. C. Polanyi, S. D. Rosner,
                                        and C. E. Young,
                                        "Energy Distribution Among Products of Exothermic Reactions. 
                                        II. Repulsive, Mixed, and Attractive Energy Release",
                                        J. Chem. Phys. 44, 1168-1184 (1966).
    ============================================================================================
    10.   ClH2_LEPS_ALTG1_1996:    single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2g3.f 
                                   ref: T. C. Allison, G. C. Lynch, D. G. Truhlar, and M. S. Gordon,
                                        "An Improved Potential Energy Surface for the H2Cl System and 
                                        Its Use for Calculations of Rate Coefficients and Kinetic Isotope 
                                        Effects",
                                        J. Phys. Chem. 100, 13575-13587 (1996).
    ============================================================================================
    11.  ClH2_LEPS_ALTG2_1996:     single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2g3v.f 
                                   ref: T. C. Allison, G. C. Lynch, D. G. Truhlar, and M. S. Gordon,
                                        "An Improved Potential Energy Surface for the H2Cl System and 
                                        Its Use for Calculations of Rate Coefficients and Kinetic Isotope 
                                        Effects",
                                        J. Phys. Chem. 100, 13575-13587 (1996).
    ============================================================================================
    12.  ClH2_LEPS_SPK_1973:       single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2spkgsw.f
                                   ref: M. J. Stern, A. Persky, and F. S. Klein,
                                        "Force field and tunneling effects in the H–H–Cl reaction system. 
                                        Determination from kinetic‐isotope‐effect",
                                        J. Phys. Chem. 58, 5697-5706 (1973).
    ============================================================================================
    13.  ClH2_LEPS_LB_1980:        single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2dim3c.f
                                   ref: I. Last, and M. Baer, 
                                        "Semi-empirical potential energy surfaces for the reactions 
                                        H + HCl → H2 + Cl and H' + HCl → H' Cl + H", 
                                        Chem. Phys. Lett. 73, 514-518 (1980).
    ============================================================================================
    14.  ClH2_LEPS_TTGI1_1985:     single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2dima.f 
                                   ref: S. C. Tucker, D. G. Truhlar, G. C. Garrett, and A. D. Isaacson,
                                        "Variational transition state theory with least‐action tunneling 
                                        calculations for the kinetic isotope effects in the Cl+H2 reaction: 
                                        Tests of extended‐LEPS, information‐theoretic, and diatomics‐in
                                        ‐molecules potential energy surfaces",
                                        J. Chem. Phys. 82, 4102-4119 (1985).
    ============================================================================================
    15.  ClH2_LEPS_TTGI2_1985:     single-state surface of ClH2, ground state,
                                   availability: potential energy
                                   functional form: LEPS
                                   corresponding surface in POTLIB: clh2dims.f 
                                   ref: S. C. Tucker, D. G. Truhlar, G. C. Garrett, and A. D. Isaacson,
                                        "Variational transition state theory with least‐action tunneling 
                                        calculations for the kinetic isotope effects in the Cl+H2 reaction: 
                                        Tests of extended‐LEPS, information‐theoretic, and diatomics‐in
                                        ‐molecules potential energy surfaces",
                                        J. Chem. Phys. 82, 4102-4119 (1985).

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_chlorine=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Cl'):
                 n_chlorine=n_chlorine+1
        if n_hydrogen!=2:
            print("number of Hydrogen atoms not equal 2 for ClH2 system")
            sys.exit()
        if n_chlorine!=1:
            print("number of Chlorine atoms not equal 1 for ClH2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)
     
    if surface=='ClH2_LEPS_GTM1_1981' or surface=='ClH2_LEPS_GTM2_1981' or surface=='ClH2_LEPS_GTM3_1981' or surface=='ClH2_LEPS_GTM4_1982':
        geom_ordered=sorted(geom, key=lambda x: ("Cl", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Cl", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    if surface=='ClH2_LEPS_PK1_1966' or surface=='ClH2_LEPS_PK2_1966':
        geom_ordered=sorted(geom, key=lambda x: ("Cl", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Cl", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    if surface=='ClH2_LEPS_KNPRY_1966' or surface=='ClH2_LEPS_STSBLTG1_1989' or surface=='ClH2_LEPS_STSBLTG2_1989':
        geom_ordered=sorted(geom, key=lambda x: ("Cl", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Cl", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    if surface=='ClH2_LEPS_ALTG1_1996' or surface=='ClH2_LEPS_ALTG2_1996' or surface=='ClH2_LEPS_SPK_1973':
        geom_ordered=sorted(geom, key=lambda x: ("Cl", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Cl", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]

    if surface=='ClH2_LEPS_LB_1980' or surface=='ClH2_LEPS_TTGI1_1985' or surface=='ClH2_LEPS_TTGI2_1985':
        geom_ordered=sorted(geom, key=lambda x: ("Cl", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Cl", "H").index(x[0]))
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

