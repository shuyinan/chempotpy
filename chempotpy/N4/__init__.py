'''
N4
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for N4:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga
    1.  N4_1A_ZV:            single-state surface of N4, singlet, ground state,      P/G
    2.  N4_BDTC_2014:        single-state surface of N4, singlet, ground state,      P/G

    Single-State Surfaces:
    1.  N4_1A_PIPNN:         single-state surface of N4, singlet, ground state,      P/G
    2.  N4_1A_PIP:           single-state surface of N4, singlet, ground state,      P/G
    3.  N4_1A_PIP_2013:      single-state surface of N4, singlet, ground state,      P/G


        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for N4:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga
    1.  N4_1A_ZV:            single-state surface of N4, singlet, ground state,
                             availability: potential energy, gradient
                             functional form: PIP
                             corresponding surface in POTLIB: N4_1A_MB-PIP-MEG3.f90
                             ref: J. Li, Z. Varga, D. G. Truhlar, and H. Guo,
                                  "Many-Body Permutationally Invariant Polynomial
                                  Neural Network Potential Energy Surface for N4",
                                  J. Chem. Theory Comput. 16, 4822-4832 (2020).
    ============================================================================================
    2.  N4_BDTC_2014:        single-state surface of N4, singlet, ground state,
                             availability: potential energy, gradient
                             functional form: PIP
                             corresponding surface in POTLIB: n4pes_L-IMLS-G2_v1.tar
                             ref: J. D. Bender, S. Doraiswamy, D. G. Truhlar, and G. V. Candler,
                                  "Potential energy surface fitting by a statistically localized, 
                                  permutationally invariant, local interpolating moving least 
                                  squares method for the many-body potential: Method and application 
                                  to N4",
                                  J. Chem. Phys. 140, 054302 (2014).


    Single-State Surfaces:
    1.  N4_1A_PIPNN:         single-state surface of N4, singlet, ground state,
                             availability: potential energy, gradient
                             functional form: PIP-NN
                             corresponding surface in POTLIB: PES_N4_MB-PIP-NN.F 
                             ref: J. Li, Z. Varga, D. G. Truhlar, and H. Guo,
                                  "Many-Body Permutationally Invariant Polynomial
                                  Neural Network Potential Energy Surface for N4",
                                  J. Chem. Theory Comput. 16, 4822-4832 (2020).
    ============================================================================================
    2.  N4_1A_PIP:           single-state surface of N4, singlet, ground state,
                             availability: potential energy, gradient
                             functional form: PIP
                             corresponding surface in POTLIB: PES_N4_singlet_umn_v3.f90
                             ref: J. D. Bender, P. Valentini, I. Nompelis, Y. Paukku, 
                                  Z. Varga, D. G. Truhlar, T. Schwartzentruber, and 
                                  G. V. Candler, 
                                  "An improved potential energy surface and multi-temperature 
                                  quasiclassical trajectory calculations of N2 + N2 dissociation 
                                  reactions"
                                  J. Chem. Phys. 143, 054304 (2015).
    ============================================================================================
    2.  N4_1A_PIP_2013:      single-state surface of N4, singlet, ground state,
                             availability: potential energy, gradient
                             functional form: PIP
                             corresponding surface in POTLIB: n4pes-gpip-v2.f
                             ref: Y. Paukku, K. R. Yang, Z, Varga, and D. G. Truhlar,
                                  "Global ab initio ground-state potential energy surface of N4",
                                  J. Chem. Phys. 139, 044309 (2013). 

        """)

    return

def check(system, surface, geom):

    n_nitrogen=0
    natoms=len(geom)

    if natoms!=4:
        print("number of atoms not equal 4 for N4 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='N'):
                n_nitrogen=n_nitrogen+1
        if n_nitrogen!=4:
            print("number of Nitrogen atoms not equal 4 for N4 system")
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
