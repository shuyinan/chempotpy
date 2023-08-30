'''
H4O2
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for H4O2:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   H4O2_PIP:         single-state surface of H4O2, ground state,      P


        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for H4O2:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   H4O2_PIP:         single-state surface of H4O2, ground state,
                           availability: potential energy
                           functional form: PIP
                           corresponding surface in POTLIB: H4O2.zip
                           ref: X. Huang, B. J. Braams, J. M. Bowman,
                                "Ab Initio Potential Energy and Dipole Moment Surfaces 
                                of (H2O)2",
                                J. Phys. Chem. A 110, 445-451 (2006).


        """)

    return

def check(system, surface, geom):

    n_oxygen=0
    n_hydrogen=0
    natoms=len(geom)

    if natoms!=6:
        check=0
        print("number of atoms not equal 6")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='O'):
                 n_oxygen=n_oxygen+1
        if n_hydrogen!=4:
            print("number of Hydrogen atoms not equal 4 for H4O2 system")
            sys.exit()
        if n_oxygen!=2:
            print("number of Oxygen atoms not equal 2 for H4O system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for H4O2_PIP: input Cartesian should in order of HHHHOO
    if surface=='H4O2_PIP':
        geom_ordered=sorted(geom, key=lambda x: ("H", "O").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "O").index(x[0]))
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
