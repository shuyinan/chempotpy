'''
N2O2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for N2O2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga
    1.  N2O2_3A_ZV:          single-state surface of N2O2, triplet, ground state,      P/G

    Single-State Surfaces:
    1.  N2O2_3A_PIP:         single-state surface of N2O2, triplet, ground state,      P/G


        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for N2O2:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    ===using unified N2, O2, and NO pairwise potentials developed by Z. Varga
    1.  N2O2_3A_ZV:          single-state surface of N2O2, triplet, ground state,
                             availability: potential energy, gradient
                             functional form: PIP
                             corresponding surface in POTLIB: N2O2_3A_MB-PIP-MEG2.f90
                             ref: Z. Varga, R. Meana-Pañeda, G. Song, Y. Paukku,
                                  and D. G. Truhlar,
                                  "Potential energy surface of triplet N2O2"
                                  J. Chem. Phys. 144, 024310 (2016).


    Single-State Surfaces:
    1.  N2O2_3A_PIP:         single-state surface of N2O2, triplet, ground state,
                             availability: potential energy, gradient
                             functional form: PIP
                             corresponding surface in POTLIB: PES_N2O2_triplet_umn_v1.f90
                             ref: Z. Varga, R. Meana-Pañeda, G. Song, Y. Paukku,
                                  and D. G. Truhlar,
                                  "Potential energy surface of triplet N2O2"
                                  J. Chem. Phys. 144, 024310 (2016).


        """)

    return

def check(system, surface, geom):

    n_oxygen=0
    n_nitrogen=0
    natoms=len(geom)

    if natoms!=4:
        print("number of atoms not equal 4 for N2O2 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='O'):
                n_oxygen=n_oxygen+1
            elif (geom[iatom][0]=='N'):
                 n_nitrogen=n_nitrogen+1
        if n_oxygen!=2:
            print("number of Oxygen atoms not equal 2 for N2O2 system")
            sys.exit()
        if n_nitrogen!=2:
            print("number of Nitrogen atoms not equal 2 for N2O2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for N2O2_3A_ZV, N2O2_3A_PIP: input Cartesian should in order of OONN
    if surface=='N2O2_3A_ZV' or surface=='N2O2_3A_PIP':
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
