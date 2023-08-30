'''
HOCO_Ar
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for HOCO_Ar:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   HOCO_Ar_PIP:       single states surface of HOCO_Ar, ground state,     P

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for HOCO_Ar:
    GEN: General
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    ============================================================================================
    1.   HOCO_Ar_PIP:       single states surface of HOCO_Ar, ground state,
                            availability: potential energy
                            functional form: PIP
                            corresponding surface in POTLIB: HOCO-Ar/HOCO-Ar.zip
                            ref: R. Conte, P. L. Houston, J. M. Bowman,
                                 "Communication: A benchmark-quality, full-dimensional 
                                 ab initio potential energy surface for Ar-HOCO",
                                 J. Chem. Phys. 140, 151101 (2014).

        """)

    return

def check(system, surface, geom):

    n_argon=0
    n_carbon=0
    n_oxygen=0
    n_hydrogen=0
    natoms=len(geom)

    if natoms!=5:
        check=0
        print("number of atoms not equal 3")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='Ar'):
                 n_argon=n_argon+1
            elif (geom[iatom][0]=='C'):
                 n_carbon=n_carbon+1
            elif (geom[iatom][0]=='O'):
                 n_oxygen=n_oxygen+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1

        if n_argon!=1:
            print("number of Argon atoms not equal 1 for HOCO_Ar system")
            sys.exit()
        if n_carbon!=1:
            print("number of Carbon atoms not equal 1 for HOCO_Ar system")
            sys.exit()
        if n_oxygen!=2:
            print("number of Oxygen atoms not equal 2 for HOCO_Ar system")
            sys.exit()
        if n_hydrogen!=1:
            print("number of Hydrogen atoms not equal 1 for HOCO_Ar system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for HOCO_Ar_PIP: input Cartesian should in order of OOCHAr
    if surface=='HOCO_Ar_PIP':
        geom_ordered=sorted(geom, key=lambda x: ("O", "C", "H", "Ar").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("O", "C", "H", "Ar").index(x[0]))
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
