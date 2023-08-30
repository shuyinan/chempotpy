'''
K2Rb2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for K2Rb2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  K2Rb2_PIPNN:         single-state surface of K2Rb2, ground state,      P



        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for K2Rb2:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  K2Rb2_PIPNN:         single-state surface of K2Rb2, ground state,
                             availability: potential energy
                             functional form: PIP-NN
                             corresponding surface in POTLIB: K2Rb2.zip
                             ref: D. Yang, J. Zuo, J. Huang, X. Hu, R. Dawes, D. Xie,
                                  and H. Guo, 
                                  "A Global Full-Dimensional Potential Energy Surface 
                                  for theK2Rb2 Complex and Its Lifetime",
                                  J. Phys. Chem. Lett. 11, 2605−2610 (2020).



        """)

    return

def check(system, surface, geom):

    n_potassium=0
    n_rubidium=0
    natoms=len(geom)

    if natoms!=4:
        print("number of atoms not equal 4 for K2Rb2 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='K'):
                n_potassium=n_potassium+1
            elif (geom[iatom][0]=='Rb'):
                 n_rubidium=n_rubidium+1
        if n_potassium!=2:
            print("number of Potassium atoms not equal 2 for K2Rb2 system")
            sys.exit()
        if n_rubidium!=2:
            print("number of Nitrogen atoms not equal 2 for K2Rb2 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for K2Rb2_PIPNN: input Cartesian should in order of KKRbRb
    if surface=='K2Rb2_PIPNN':
        geom_ordered=sorted(geom, key=lambda x: ("K", "Rb").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("K", "Rb").index(x[0]))
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
