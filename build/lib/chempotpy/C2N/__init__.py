'''
C2N
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for C2N:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  C2N_PIPNN_Ap:         single-state surface of C2N, A' symmetry ground state,   P
    2.  C2N_PIPNN_App:        single-state surface of C2N, A'' symmetry ground state,  P

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule: 
    The ordering will be CCN     

        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for C2N:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  C2N_PIPNN_Ap:         single-state surface of C2N, A' symmetry ground state,   P
                              availability: potential energy
                              functional form: PIP-NN
                              corresponding surface in POTLIB: c2n_2a_prime.f 
                              ref: J. Zuo, D. Zhang, D. G. Truhlar, H. Guo, 
                                   "Global Potential Energy Surfaces by Compressed-State 
                                   Multistate Pair-Density Functional Theory: The Lowest 
                                   Doublet States Responsible for the N(4Su) + C2(a 3Πu) 
                                   → CN(X 2Σ+) + C(3Pg) Reaction"
                                   J. Chem. Theory Comput. 18, 7121-7131 (2022).
    2.  C2N_PIPNN_App:        single-state surface of C2N, A'' symmetry ground state,  P
                              availability: potential energy
                              functional form: PIP-NN
                              corresponding surface in POTLIB: c2n_2a_prime.f 
                              ref: J. Zuo, D. Zhang, D. G. Truhlar, H. Guo, 
                                   "Global Potential Energy Surfaces by Compressed-State 
                                   Multistate Pair-Density Functional Theory: The Lowest 
                                   Doublet States Responsible for the N(4Su) + C2(a 3Πu) 
                                   → CN(X 2Σ+) + C(3Pg) Reaction"
                                   J. Chem. Theory Comput. 18, 7121-7131 (2022).

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule:
    The ordering will be CCN  

        """)

    return

def check(system, surface, geom):

    n_nitrogen=0
    n_carbon=0
    natoms=len(geom)

    if natoms!=3:
        print("number of atoms not equal 3 for C2N system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='N'):
                 n_nitrogen=n_nitrogen+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
        if n_carbon!=2:
            print("number of Carbon atoms not equal 2 for C2N system")
            sys.exit()
        if n_nitrogen!=1:
            print("number of Hydrogen atoms not equal 1 for C2N system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    if surface=='C2N_PIPNN_Ap' or surface=='C2N_PIPNN_App': 
        geom_ordered=sorted(geom, key=lambda x: ("C", "N").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("C", "N").index(x[0]))
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
