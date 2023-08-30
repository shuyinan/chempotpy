'''
NO2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for NO2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga
    1.  NO2_2Ap_ZV:         single-state surface of NO2, doublet, A' symmetry ground state,      P/G
    2.  NO2_4Ap_ZV:         single-state surface of NO2, quartet, A' symmetry ground state,      P/G
    3.  NO2_6Ap_ZV:         single-state surface of NO2, sextet, A' symmetry ground state,       P/G


    Single-State Surfaces:
    1.  NO2_2Ap_PIPNN:      single-state surface of NO2, doublet, A' symmetry ground state,      P/G
    2.  NO2_4Ap_PIPNN:      single-state surface of NO2, quartet, A' symmetry ground state,      P/G
    3.  NO2_6Ap_PIPNN:      single-state surface of NO2, sextet, A' symmetry ground state,       P/G


        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for NO2:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga
    1.  NO2_2Ap_ZV:          single-state surface of NO2, doublet, A' symmetry ground state,
                             availability: potential energy, gradient
                             functional form: PIP
                             corresponding surface in POTLIB: NO2_2Ap_MB-PIP-MEG.f90
                             ref: Z. Varga, Y. Liu, J. Li, Y. Paukku, H. Guo, and D. G. Truhlar
                             "Potential energy surfaces for high-energy N + O2 collisions",
                             J. Chem. Phys. 154, 084304 (2021). 
    ============================================================================================
    2.  NO2_4Ap_ZV:          single-state surface of NO2, quartet, A' symmetry ground state,
                             availability: potential energy, gradient
                             functional form: PIP
                             corresponding surface in POTLIB: NO2_4Ap_MB-PIP-MEG.f90 
                             ref: Z. Varga, Y. Liu, J. Li, Y. Paukku, H. Guo, and D. G. Truhlar
                             "Potential energy surfaces for high-energy N + O2 collisions",
                             J. Chem. Phys. 154, 084304 (2021). 
     ============================================================================================
    3.  NO2_4Ap_ZV:          single-state surface of NO2, sextet, A' symmetry ground state,
                             availability: potential energy, gradient
                             functional form: PIP
                             corresponding surface in POTLIB: NO2_6Ap_MB-PIP-MEG.f90 
                             ref: Z. Varga, Y. Liu, J. Li, Y. Paukku, H. Guo, and D. G. Truhlar
                             "Potential energy surfaces for high-energy N + O2 collisions",
                             J. Chem. Phys. 154, 084304 (2021). 

    Single-State Surfaces:
    1.  NO2_2Ap_PIPNN:       single-state surface of NO2, doublet, A' symmetry ground state,
                             availability: potential energy, gradient
                             functional form: PIP-NN
                             corresponding surface in POTLIB: NO2_2Ap_PIP-NN.zip 
                             ref: Z. Varga, Y. Liu, J. Li, Y. Paukku, H. Guo, and D. G. Truhlar
                             "Potential energy surfaces for high-energy N + O2 collisions",
                             J. Chem. Phys. 154, 084304 (2021). 
    ============================================================================================
    2.  NO2_4Ap_PIPNN:       single-state surface of NO2, quartet, A' symmetry ground state,
                             availability: potential energy, gradient
                             functional form: PIP-NN
                             corresponding surface in POTLIB: NO2_4Ap_PIP-NN.zip
                             ref: Z. Varga, Y. Liu, J. Li, Y. Paukku, H. Guo, and D. G. Truhlar
                             "Potential energy surfaces for high-energy N + O2 collisions",
                             J. Chem. Phys. 154, 084304 (2021). 
     ============================================================================================
    3.  NO2_4Ap_PIPNN:       single-state surface of NO2, sextet, A' symmetry ground state,
                             availability: potential energy, gradient
                             functional form: PIP-NN
                             corresponding surface in POTLIB: NO2_6Ap_PIP-NN.zip
                             ref: Z. Varga, Y. Liu, J. Li, Y. Paukku, H. Guo, and D. G. Truhlar
                             "Potential energy surfaces for high-energy N + O2 collisions",
                             J. Chem. Phys. 154, 084304 (2021). 


        """)


    return

def check(system, surface, geom):

    n_oxygen=0
    n_nitrogen=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3 for NO2 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                n_oxygen=n_oxygen+1
            elif (geom[iatom][0]=='N'):
                 n_nitrogen=n_nitrogen+1
        if n_oxygen!=2:
            print("number of Oxygen atoms not equal 2 for NO2 system")
            sys.exit()
        if n_nitrogen!=1:
            print("number of Nitrogen atoms not equal 1 for NO2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for NO2_2Ap_ZV, NO2_4Ap_ZV, NO2_6Ap_ZV: input Cartesian should in order of OON
    if surface=='NO2_2Ap_ZV' or surface=='NO2_4Ap_ZV' or surface=='NO2_6Ap_ZV' or surface=='NO2_2Ap_PIPNN' or surface=='NO2_4Ap_PIPNN' or surface=='NO2_6Ap_PIPNN': 
        geom_ordered=sorted(geom, key=lambda x: (x[0], x[1]),reverse=True)
        rank_ordered=sorted(rank, key=lambda x: (x[0], x[1]),reverse=True)
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
