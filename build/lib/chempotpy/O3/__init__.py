'''
O3
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for O3:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    U: diabatic potential energy matrix 
    UG: gradient of diabatic potential energy matrix
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Recommended Surfaces for Use:
    M1.  O3_14_3Ap_2023:   multi-state surface of O3, triplet, A' symmetry 1-14 state,         P/G/D
    M2.  O3_6_5Ap_2023:    multi-state surface of O3, quintet, A' symmetry 1-6 state,          P/G/D
    Z1.  O3_1_1Ap_ZV:      single-state surface of O3, singlet, A' symmetry ground state,      P/G
    Z2.  O3_1_1App_ZV:     single-state surface of O3, singlet, A'' symmetry ground state,     P/G
    Z3.  O3_1_3Ap_ZV:      single-state surface of O3, triplet, A' symmetry ground state,      P/G
    Z4.  O3_1_3App_ZV:     single-state surface of O3, triplet, A'' symmetry ground state,     P/G
    Z5.  O3_1_5Ap_ZV:      single-state surface of O3, quintet, A' symmetry ground state,      P/G
    Z6.  O3_1_5App_ZV:     single-state surface of O3, quintet, A'' symmetry ground state,     P/G
    Z7.  O3_2_1Ap_ZV:      single-state surface of O3, singlet, A' symmetry 1st excited state, P/G
    Z8.  O3_2_3Ap_ZV:      single-state surface of O3, triplet, A' symmetry 1st excited state, P/G
    Z9.  O3_2_5Ap_ZV:      single-state surface of O3, quintet, A' symmetry 1st excited state, P/G


    Complete List of Surfaces:
    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga 
    Z1.  O3_1_1Ap_ZV:      single-state surface of O3, singlet, A' symmetry ground state,      P/G
    Z2.  O3_1_1App_ZV:     single-state surface of O3, singlet, A'' symmetry ground state,     P/G
    Z3.  O3_1_3Ap_ZV:      single-state surface of O3, triplet, A' symmetry ground state,      P/G
    Z4.  O3_1_3App_ZV:     single-state surface of O3, triplet, A'' symmetry ground state,     P/G
    Z5.  O3_1_5Ap_ZV:      single-state surface of O3, quintet, A' symmetry ground state,      P/G
    Z6.  O3_1_5App_ZV:     single-state surface of O3, quintet, A'' symmetry ground state,     P/G
    Z7.  O3_2_1Ap_ZV:      single-state surface of O3, singlet, A' symmetry 1st excited state, P/G
    Z8.  O3_2_3Ap_ZV:      single-state surface of O3, triplet, A' symmetry 1st excited state, P/G
    Z9.  O3_2_5Ap_ZV:      single-state surface of O3, quintet, A' symmetry 1st excited state, P/G


    Single-State Surfaces:
    S1.  O3_1_1Ap:         single-state surface of O3, singlet, A' symmetry ground state,      P/G
    S2.  O3_1_1App:        single-state surface of O3, singlet, A'' symmetry ground state,     P/G
    S3.  O3_1_3Ap:         single-state surface of O3, triplet, A' symmetry ground state,      P/G
    S4.  O3_1_3App:        single-state surface of O3, triplet, A'' symmetry ground state,     P/G
    S5.  O3_1_5Ap:         single-state surface of O3, quintet, A' symmetry ground state,      P/G
    S6.  O3_1_5App:        single-state surface of O3, quintet, A'' symmetry ground state,     P/G
    S7.  O3_2_1Ap:         single-state surface of O3, singlet, A' symmetry 1st excited state, P/G
    S8.  O3_2_3Ap:         single-state surface of O3, triplet, A' symmetry 1st excited state, P/G
    S9.  O3_2_5Ap:         single-state surface of O3, quintet, A' symmetry 1st excited state, P/G
    S10. O3_MBP_2017:      single-state surface of O3, ground state,                           P


    Multi-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga 
    MZ1.  O3_14_3Ap_2022:       multi-state surface of O3, triplet, A' symmetry 1-14 state,    P/G/D
    MZ2.  O3_14_3Ap_2022_DPEM:  multi-state surface of O3, triplet, A' symmetry 1-14 state,    U/UG


    Multi-State Surfaces:
    ===using a new O2 pairwise potentials improved by F. B. Ahker
    M1.  O3_14_3Ap_2023:        multi-state surface of O3, triplet, A' symmetry 1-14 state,    P/G/D
    M2.  O3_6_5Ap_2023:         multi-state surface of O3, quintet, A' symmetry 1-6 state,     P/G/D
    M3.  O3_14_3Ap_2023_DPEM:   multi-state surface of O3, triplet, A' symmetry 1-14 state,    U/UG
    M4.  O3_6_5Ap_2023_DPEM:    multi-state surface of O3, quintet, A' symmetry 1-6 state,     U/UG


        """)

    return


def intro_detail():
 
    print("""
    DETAILED List of potential energy surfaces for O3:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga 
    1.  O3_1_1Ap_ZV:      single-state surface of O3, singlet, A' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: PES_O3_1_1Ap_umn_v1.f90 
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    2.  O3_1_1App_ZV:     single-state surface of O3, singlet, A'' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: PES_O3_1_1App_umn_v1.f90
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    3.  O3_1_1App_ZV:     single-state surface of O3, triplet, A' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: PES_O3_1_3Ap_umn_v1.f90
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    4.  O3_1_3App_ZV:     single-state surface of O3, triplet, A'' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: PES_O3_1_3App_umn_v1.f90
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    5.  O3_1_5Ap_ZV:      single-state surface of O3, quintet, A' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: PES_O3_1_5Ap_umn_v1.f90
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    6.  O3_1_5App_ZV:     single-state surface of O3, quintet, A'' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: PES_O3_1_5App_umn_v1.f90 
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    7.  O3_2_1Ap_ZV:      single-state surface of O3, singlet, A' symmetry 1st excited state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: PES_O3_2_1Ap_umn_v1.f90
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    8.  O3_2_3Ap_ZV:      single-state surface of O3, triplet, A' symmetry 1st excited state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: PES_O3_2_3Ap_umn_v1.f90
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    9.  O3_2_5Ap_ZV:      single-state surface of O3, quintet, A' symmetry 1st excited state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: PES_O3_2_5Ap_umn_v1.f90
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).



    Single-State Surfaces:
    1.  O3_1_1Ap:         single-state surface of O3, singlet, A' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: O3_1_1Ap.f
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    2.  O3_1_1App:        single-state surface of O3, singlet, A'' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: O3_1_1App.f
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    3.  O3_1_1App:        single-state surface of O3, triplet, A' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: O3_1_3Ap.f 
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    4.  O3_1_3App:        single-state surface of O3, triplet, A'' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: O3_1_3App.f
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    5.  O3_1_5Ap:         single-state surface of O3, quintet, A' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: O3_1_5Ap.f  
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    6.  O3_1_5App:        single-state surface of O3, quintet, A'' symmetry ground state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: O3_1_5App.f
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    7.  O3_2_1Ap:         single-state surface of O3, singlet, A' symmetry 1st excited state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: O3_2_1Ap.f 
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    8.  O3_2_3Ap:         single-state surface of O3, triplet, A' symmetry 1st excited state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: O3_2_3Ap.f
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    9.  O3_2_5Ap:         single-state surface of O3, quintet, A' symmetry 1st excited state,
                          availability: potential energy, gradient 
                          functional form: PIP
                          corresponding surface in POTLIB: O3_2_5Ap.f 
                          ref: Z. Varga, Y. Paukku, and D. G. Truhlar,
                               "Potential energy surfaces for O + O2 collisions",
                               J. Chem. Phys. 147, 154312/1-17 (2017).
    ==================================================================================== 
    10. O3_MBP_2017:      single-state surface of O3, ground electronic state; but spin and
                          spatial symmetry information is unknow. 
                          availability: potential energy
                          functional form: polynomial
                          corresponding surface in POTLIB: o3pes.f90 
                          ref: T. K. Mankodi, U. V. Bhandarkar, B. P. Puranik,
                               "Dissociation cross sections for N2+N → 3N and O2+O → 3O using
                               the QCT method",
                               J. Chem. Phys. 146, 204307 (2017).


    Multi-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga
    1.  O3_14_3Ap_2022:   multi-state surface of O3, triplet, A' symmetry 1-14 state,
                          availability: potential energy, gradient, nonadiabatic coupling
                          functional form: PR-DDNN
                          corresponding surface in POTLIB: O3_14_3Ap_2022.f90
                          ref: Z. Varga, Y. Shu, J. Ning, D. G. Truhlar, 
                               "Diabatic potential energy surfaces and semiclassical multi-state
                               dynamics for fourteen coupled 3A' states of O3",
                               Electron. Struct. 4, 047002 (2022).



    Multi-State Surfaces:
    ===using a new O2 pairwise potentials improved by F. B. Ahker
    1.  O3_14_3Ap_2023:   multi-state surface of O3, triplet, A' symmetry 1-14 state,
                          availability: potential energy, gradient, nonadiabatic coupling
                          functional form: PM-DDNN 
                          corresponding surface in POTLIB: O3_14_3Ap_2023.f90 
                          ref: 

    ==================================================================================== 
    2.  O3_6_5Ap_2023:   multi-state surface of O3, quintet, A' symmetry 1-6 state, 
                         availability: potential energy, gradient, nonadiabatic coupling
                         functional form: PM-DDNN 
                          corresponding surface in POTLIB: O3_6_5Ap_2023.f90
                         ref:  



        """)


    return


def check(system, surface, geom):

    n_oxygen=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3 for O3 system")
        sys.exit() 
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                n_oxygen=n_oxygen+1
        if n_oxygen!=3:
            print("number of Oxygen atoms not equal 3 for O3 system")
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
