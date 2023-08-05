'''
H3p
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for H3p:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   H3p_UNO_2001      single-state surface of H3p, ground state,  P/G
    2.   H3p_KBNN_2002     single-state surface of H3p, ground state,  P/G

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for H3p:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   H3p_UNO_2001      multi-state surface of H3p, ground state,
                           availability: potential energy, gradient
                           functional form: GEN
                           corresponding surface in POTLIB: h3plusUNOd.f 
                           ref: V. G. Ushakov, K. Nobusada, V. I. Osherov, 
                                "Electronically nonadiabatic transitions in a collinear H2 + H‘ 
                                system: Quantum mechanical understanding and comparison with a 
                                trajectory surface hopping method",
                                Phys. Chem. Chem. Phys. 3, 63-69 (2000).
    ============================================================================================
    2.   H3p_KBNN_2002     multi-state surface of H3p, ground state,
                           availability: potential energy, gradient
                           functional form: GEN
                           corresponding surface in POTLIB: h3plusKBNNd.f 
                           ref: H. Kamisaka, W. Bian, K. Nobusada, and H. Nakamura,
                                "Accurate quantum dynamics of electronically nonadiabatic chemical 
                                reactions in the DH+2 system", 
                                J. Chem. Phys. 116, 654-665 (2002).


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
            print("number of Hydrogen atoms not equal 3 for H3p system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='H3p_UNO_2001' or surface=='H3p_KBNN_2002':
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
