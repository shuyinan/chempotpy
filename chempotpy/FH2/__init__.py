'''
FH2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for FH2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   FH2_GEN_5SEC_1991:       single-state surface of FH2, ground state,      P/G
    2.   FH2_GEN_5SECCOG1_1993:   single-state surface of FH2, ground state,      P/G
    3.   FH2_GEN_5SECCOG2_1993:   single-state surface of FH2, ground state,      P/G
    4.   FH2_GEN_5SECCOG3_1993:   single-state surface of FH2, ground state,      P/G
    5.   FH2_GEN_5SECG1_1993:     single-state surface of FH2, ground state,      P/G
    6.   FH2_GEN_5SECG2_1993:     single-state surface of FH2, ground state,      P/G
    7.   FH2_GEN_5SECG3_1993:     single-state surface of FH2, ground state,      P/G
    8.   FH2_GEN_5SECS1_1993:     single-state surface of FH2, ground state,      P/G
    9.   FH2_GEN_5SECS2_1993:     single-state surface of FH2, ground state,      P/G
    10.  FH2_GEN_5SECS3_1993:     single-state surface of FH2, ground state,      P/G
    11.  FH2_GEN_5SECW_1991:      single-state surface of FH2, ground state,      P/G
    12.  FH2_GEN_6SEC_1993:       single-state surface of FH2, ground state,      P/G
    13.  FH2_LEPS_M5_1991:        single-state surface of FH2, ground state,      P/G
    14.  FH2_LEPS_M5SK_1966:      single-state surface of FH2, ground state,      P/G
    15.  FH2_LEPS_NO1_1984:       single-state surface of FH2, ground state,      P/G
    16.  FH2_LEPS_NO2_1984:       single-state surface of FH2, ground state,      P/G
    17.  FH2_LEPS_NO2a_1984:      single-state surface of FH2, ground state,      P/G
    18.  FH2_LEPS_NO3_1984:       single-state surface of FH2, ground state,      P/G
    19.  FH2_LEPS_NO4_1984:       single-state surface of FH2, ground state,      P/G
    20.  FH2_LEPS_NO5_1985:       single-state surface of FH2, ground state,      P/G
    21.  FH2_LEPS_NO5a_1985:      single-state surface of FH2, ground state,      P/G
    22.  FH2_LEPS_NO5b_1985:      single-state surface of FH2, ground state,      P/G
    23.  FH2_LEPS_NO5z_1985:      single-state surface of FH2, ground state,      P/G
    24.  FH2_LEPS_SK_1980:        single-state surface of FH2, ground state,      P/G 

        """)

    return


def intro_detail():
 
    print("""
    DETAILED List of potential energy surfaces for F2H:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   FH2_GEN_5SEC_1991:       single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh25sec2001.f 
                                  ref: G. C. Lynch, R. Steckler, D. W. Schwenke,
                                       A. J. C. Varandas, D. G. Truhlar, and B. C. Garrett,
                                       "Use of scaled external correlation, a double many‐
                                       body expansion, and variational transition state theory
                                       to calibrate a potential energy surface for FH2",
                                       J. Chem. Phys. 94, 7136-7419 (1991). 
    ============================================================================================
    2.   FH2_GEN_5SECCOG1_1993:   single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh25seccog1-2001.f 
                                  ref: S. L. Mielke, G. C. Lynch, D. G. Truhlar, and D. W. Schwenke,
                                       "A more accurate potential energy surface and quantum mechanical 
                                       cross section calculations for the F+ H2 reaction", 
                                       Chem. Phys. Lett. 213, 10-16 (1993).
    ============================================================================================
    3.   FH2_GEN_5SECCOG2_1993:   single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh25seccog2-2001.f 
                                  ref: S. L. Mielke, G. C. Lynch, D. G. Truhlar, and D. W. Schwenke,
                                       "A more accurate potential energy surface and quantum mechanical 
                                       cross section calculations for the F+ H2 reaction", 
                                       Chem. Phys. Lett. 213, 10-16 (1993).
    ============================================================================================
    4.   FH2_GEN_5SECCOG3_1993:   single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh25seccog3-2001.f 
                                  ref: S. L. Mielke, G. C. Lynch, D. G. Truhlar, and D. W. Schwenke,
                                       "A more accurate potential energy surface and quantum mechanical 
                                       cross section calculations for the F+ H2 reaction", 
                                       Chem. Phys. Lett. 213, 10-16 (1993).
    ============================================================================================
    5.   FH2_GEN_5SECG1_1993:     single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh25secgm1-2001.f  
                                  ref: S. L. Mielke, G. C. Lynch, D. G. Truhlar, and D. W. Schwenke,
                                       "A more accurate potential energy surface and quantum mechanical 
                                       cross section calculations for the F+ H2 reaction", 
                                       Chem. Phys. Lett. 213, 10-16 (1993).
    ============================================================================================
    6.   FH2_GEN_5SECG2_1993:     single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh25secgm2-2001.f  
                                  ref: S. L. Mielke, G. C. Lynch, D. G. Truhlar, and D. W. Schwenke,
                                       "A more accurate potential energy surface and quantum mechanical 
                                       cross section calculations for the F+ H2 reaction", 
                                       Chem. Phys. Lett. 213, 10-16 (1993).
    ============================================================================================
    7.   FH2_GEN_5SECG3_1993:     single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh25secgm3-2001.f 
                                  ref: S. L. Mielke, G. C. Lynch, D. G. Truhlar, and D. W. Schwenke,
                                       "A more accurate potential energy surface and quantum mechanical 
                                       cross section calculations for the F+ H2 reaction", 
                                       Chem. Phys. Lett. 213, 10-16 (1993).
    ============================================================================================
   8.   FH2_GEN_5SECS1_1993:      single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh25secs1-2001.f
                                  ref: S. L. Mielke, G. C. Lynch, D. G. Truhlar, and D. W. Schwenke,
                                       "A more accurate potential energy surface and quantum mechanical 
                                       cross section calculations for the F+ H2 reaction", 
                                       Chem. Phys. Lett. 213, 10-16 (1993).
    ============================================================================================
    9.   FH2_GEN_5SECS2_1993:     single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh25secs2-2001.f
                                  ref: S. L. Mielke, G. C. Lynch, D. G. Truhlar, and D. W. Schwenke,
                                       "A more accurate potential energy surface and quantum mechanical 
                                       cross section calculations for the F+ H2 reaction", 
                                       Chem. Phys. Lett. 213, 10-16 (1993).
    ============================================================================================
    10.  FH2_GEN_5SECS3_1993:     single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh25secs3-2001.f
                                  ref: S. L. Mielke, G. C. Lynch, D. G. Truhlar, and D. W. Schwenke,
                                       "A more accurate potential energy surface and quantum mechanical 
                                       cross section calculations for the F+ H2 reaction", 
                                       Chem. Phys. Lett. 213, 10-16 (1993).
    ============================================================================================
    11.  FH2_GEN_5SECW_1991:      single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh25secw2001.f 
                                  ref: G. C. Lynch, P. Halvick, M. Zhao, D. G. Truhlar, C.-h. Yu,
                                       D. J. Kouri, and D. W. Schwenke,
                                       "Converged three‐dimensional quantum mechanical reaction 
                                       probabilities for the F+H2 reaction on a potential energy 
                                       surface with realistic entrance and exit channels and comparisons 
                                       to results for three other surfaces",
                                       J. Chem. Phys. 94, 7150-7158 (1991). 
    ============================================================================================
    12.  FH2_GEN_6SEC_1993:       single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: GEN
                                  corresponding surface in POTLIB: fh26sec2001.f 
                                  ref: S. L. Mielke, G. C. Lynch, D. G. Truhlar, and D. W. Schwenke,
                                       "A more accurate potential energy surface and quantum mechanical 
                                       cross section calculations for the F+ H2 reaction", 
                                       Chem. Phys. Lett. 213, 10-16 (1993).
    ============================================================================================
    13.  FH2_LEPS_M5_1991:        single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2m5-2001.f 
                                  ref: J. T. Muckerman,
                                       Theor. Chem. Advan. Perspectives 6A, 1 (1981)
    ============================================================================================
    14.  FH2_LEPS_M5SK_1966:      single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2m5sk2001.f 
                                  ref: P. J. Kuntz, E. M. Nemth, J. C. Polanyi, S. D. Rosner, 
                                       and C. E. Young,
                                       "Energy Distribution Among Products of Exothermic Reactions. 
                                       II. Repulsive, Mixed, and Attractive Energy Release",
                                       J. Chem. Phys. 44, 1168-1184 (1966).
    ============================================================================================
    15.  FH2_LEPS_NO1_1984:       single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2n1-2001.f 
                                  ref: D. G. Truhlar, B. C. Garrett, and N. C. Blais,
                                       "Two new potential energy surfaces for the F+H2 reaction",
                                       J. Chem. Phys. 80, 232-240 (1984).
    ============================================================================================
    16.  FH2_LEPS_NO2_1984:       single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2n2-2001.f 
                                  ref: D. G. Truhlar, B. C. Garrett, and N. C. Blais,
                                       "Two new potential energy surfaces for the F+H2 reaction",
                                       J. Chem. Phys. 80, 232-240 (1984).
    ============================================================================================
    17.  FH2_LEPS_NO2a_1984:      single-state surface of FH2, ground state, 
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2n2a-2001.f 
                                  ref: D. G. Truhlar, B. C. Garrett, and N. C. Blais,
                                       "Two new potential energy surfaces for the F+H2 reaction",
                                       J. Chem. Phys. 80, 232-240 (1984).
    ============================================================================================
    18.  FH2_LEPS_NO3_1984:       single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2n3-2001.f 
                                  ref: D. G. Truhlar, B. C. Garrett, and N. C. Blais,
                                       "Two new potential energy surfaces for the F+H2 reaction",
                                       J. Chem. Phys. 80, 232-240 (1984).
    ============================================================================================
    19.  FH2_LEPS_NO4_1984:       single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2n4-2001.f 
                                  ref: D. G. Truhlar, B. C. Garrett, and N. C. Blais,
                                       "Two new potential energy surfaces for the F+H2 reaction",
                                       J. Chem. Phys. 80, 232-240 (1984).
    ============================================================================================
    20.  FH2_LEPS_NO5_1985:       single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2n5-2001.f 
                                  ref: F. B. Brown, R. Steckler, D. W. Schwenke, D. G. Truhlar, 
                                       and B. C. Garrett, 
                                       "An improved potential energy surface for F + H2 → HF + H and 
                                       H + H′F → HF + H′", 
                                       J. Chem. Phys. 82, 188-201 (1985).
    ============================================================================================
    21.  FH2_LEPS_NO5a_1985:      single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2n5a-2001.f 
                                  ref: F. B. Brown, R. Steckler, D. W. Schwenke, D. G. Truhlar, 
                                       and B. C. Garrett, 
                                       "An improved potential energy surface for F + H2 → HF + H and 
                                       H + H′F → HF + H′", 
                                       J. Chem. Phys. 82, 188-201 (1985).
    ============================================================================================
    22.  FH2_LEPS_NO5b_1985:      single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2n5b-2001.f 
                                  ref: F. B. Brown, R. Steckler, D. W. Schwenke, D. G. Truhlar, 
                                       and B. C. Garrett, 
                                       "An improved potential energy surface for F + H2 → HF + H and 
                                       H + H′F → HF + H′", 
                                       J. Chem. Phys. 82, 188-201 (1985).
    ============================================================================================
    23.  FH2_LEPS_NO5z_1985:      single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2n5z-2001.f 
                                  ref: F. B. Brown, R. Steckler, D. W. Schwenke, D. G. Truhlar, 
                                       and B. C. Garrett, 
                                       "An improved potential energy surface for F + H2 → HF + H and 
                                       H + H′F → HF + H′", 
                                       J. Chem. Phys. 82, 188-201 (1985).
    ============================================================================================
    24.  FH2_LEPS_SK_1980:        single-state surface of FH2, ground state,
                                  availability: potential energy, gradient
                                  functional form: LEPS
                                  corresponding surface in POTLIB: fh2sk2001.f
                                  ref:  G. C. Schatz, and A. Kuppermann,
                                       "Vibrational deactivation on chemically reactive potential 
                                       surfaces: An exact quantum study of a low barrier collinear 
                                       model of H + FH, D + FD, H + FD and D + FH",
                                       J. Chem. Phys. 72, 2737-2743 (1980).  


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_fluorine=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='F'):
                 n_fluorine=n_fluorine+1
        if n_hydrogen!=2:
            print("number of Hydrogen atoms not equal 2 for FH2 system")
            sys.exit()
        if n_fluorine!=1:
            print("number of Fluorine atoms not equal 1 for FH2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='FH2_GEN_5SEC_1991' or surface=='FH2_GEN_5SECCOG1_1993' or surface=='FH2_GEN_5SECCOG2_1993' or surface=='FH2_GEN_5SECCOG3_1993' or surface=='FH2_GEN_5SECG1_1993' or surface=='FH2_GEN_5SECG2_1993' or surface=='FH2_GEN_5SECG3_1993' or surface=='FH2_GEN_5SECS1_1993' or surface=='FH2_GEN_5SECS2_1993' or surface=='FH2_GEN_5SECS3_1993' or surface=='FH2_GEN_5SECW_1991' or surface=='FH2_GEN_6SEC_1993' or surface=='FH2_LEPS_M5_1991' or surface=='FH2_LEPS_M5SK_1966' or surface=='FH2_LEPS_NO1_1984' or surface=='FH2_LEPS_NO2_1984' or surface=='FH2_LEPS_NO2a_1984' or surface=='FH2_LEPS_NO3_1984' or surface=='FH2_LEPS_NO4_1984' or surface=='FH2_LEPS_NO5_1985' or surface=='FH2_LEPS_NO5a_1985' or surface=='FH2_LEPS_NO5b_1985' or surface=='FH2_LEPS_NO5z_1985' or surface=='FH2_LEPS_SK_1980': 
        geom_ordered=sorted(geom, key=lambda x: ("F", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("F", "H").index(x[0]))
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
