'''
O4
'''

import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for O4:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga 
    1.  O4_singlet_ZV:       single-state surface of O4, singlet, ground state,      P/G
    2.  O4_triplet_ZV:       single-state surface of O4, triplet, ground state,      P/G
    3.  O4_quintet_ZV:       single-state surface of O4, quintet, ground state,      P/G

    Single-State Surfaces:
    1.  O4_singlet:          single-state surface of O4, singlet, ground state,      P/G
    2.  O4_triplet_v2:       single-state surface of O4, triplet, ground state,      P/G
    3.  O4_quintet:          single-state surface of O4, quintet, ground state,      P/G
    4.  O4_MBP_singlet_2018: single-state surface of O4, singlet, ground state,      P

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for O4:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network

    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga 
    1.  O4_singlet_ZV:       single-state surface of O4, singlet, ground state,
                             availability: potential energy, gradient 
                             functional form: PIP
                             corresponding surface in POTLIB: PES_O4_singlet_umn_v1.f90
                             special emphsize on: O2 + O2 → O2 + O + O reaction 
                             ref: Y. Paukku, K. R. Yang, Z. Varga, G. Song, J. D. Bender,
                                  and D. G. Truhlar,
                                  "Potential energy surfaces of quintet and singlet O4",
                                  J. Chem. Phys. 147, 034301 (2017). 
    ==================================================================================== 
    2.  O4_triplet_ZV:       single-state surface of O4, triplet, ground state,
                             availability: potential energy, gradient 
                             functional form: PIP
                             corresponding surface in POTLIB: PES_O4_triplet_umn_v2.f90
                             special emphsize on: O2 + O2 → O2 + O + O reaction 
                             ref: Y. Paukku, Z. Varga, and D. G. Truhlar,
                                  "Potential energy surfaces of triplet O4", 
                                  J. Chem. Phys. 148, 124314 (2018). 
    ==================================================================================== 
    3.  O4_quintet_ZV:       single-state surface of O4, quintet, ground state,
                             availability: potential energy, gradient 
                             functional form: PIP
                             corresponding surface in POTLIB: PES_O4_quintet_umn_v1.f90 
                             special emphsize on: O2 + O2 → O2 + O + O reaction 
                             ref: Y. Paukku, K. R. Yang, Z. Varga, G. Song, J. D. Bender,
                                  and D. G. Truhlar,
                                  "Potential energy surfaces of quintet and singlet O4",
                                  J. Chem. Phys. 147, 034301 (2017). 

    Single-State Surfaces:
    1.  O4_singlet:          single-state surface of O4, singlet, ground state,
                             availability: potential energy, gradient 
                             functional form: PIP
                             corresponding surface in POTLIB: O4_singlet.f
                             special emphsize on: O2 + O2 → O2 + O + O reaction 
                             ref: Y. Paukku, K. R. Yang, Z. Varga, G. Song, J. D. Bender,
                                  and D. G. Truhlar,
                                  "Potential energy surfaces of quintet and singlet O4",
                                  J. Chem. Phys. 147, 034301 (2017). 
    ==================================================================================== 
    2.  O4_triplet_v2:       single-state surface of O4, triplet, ground state,
                             availability: potential energy, gradient 
                             functional form: PIP
                             corresponding surface in POTLIB: O4_triplet_v2.f  
                             special emphsize on: O2 + O2 → O2 + O + O reaction 
                             ref: Y. Paukku, Z. Varga, and D. G. Truhlar,
                                  "Potential energy surfaces of triplet O4", 
                                  J. Chem. Phys. 148, 124314 (2018). 
    ==================================================================================== 
    3.  O4_quintet:          single-state surface of O4, quintet, ground state,
                             availability: potential energy, gradient 
                             functional form: PIP
                             corresponding surface in POTLIB: O4_quintet.f
                             special emphsize on: O2 + O2 → O2 + O + O reaction 
                             ref: Y. Paukku, K. R. Yang, Z. Varga, G. Song, J. D. Bender,
                                  and D. G. Truhlar,
                                  "Potential energy surfaces of quintet and singlet O4",
                                  J. Chem. Phys. 147, 034301 (2017). 
    ====================================================================================
    4.  O4_MBP_singlet_2018: single-state surface of O4, singlet, ground state,
                             availability: potential energy
                             functional form: PIP
                             corresponding surface in POTLIB: o4pes.f90 
                             special emphsize on: O2 + O2 → O2 + O + O reaction and 
                                                  O2 + O2 → 4O reaction 
                             ref: T. J. Mankodi, U. V. Bhandarkar, and B. P. Puranik,
                                  "Global potential energy surface of ground state signlet 
                                  spin O4", 
                                  J. Chem. Phys. 148, 074305 (2018).

        """)


    return


def check(system, surface, geom):

    n_oxygen=0
    natoms=len(geom)

    if natoms!=4:
        print("number of atoms not equal 4 for O4 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                n_oxygen=n_oxygen+1
        if n_oxygen!=4:
            print("number of Oxygen atoms not equal 4 for O4 system")
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
