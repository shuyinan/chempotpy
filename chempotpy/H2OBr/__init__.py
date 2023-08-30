'''
H2OBr
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for H2OBr:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   H2OBr_PIP:         single-state surface of H2OBr, ground state,      P


        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for H2OBr:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   H2OBr_PIP:         single-state surface of H2OBr, ground state,
                            availability: potential energy
                            functional form: PIP
                            corresponding surface in POTLIB: H2OBr.tar.gz
                            special emphsize on: OH + HBr → Br + H2O reaction 
                            ref: A. G. S. de Oliveira-Filho, F. R. Ornellas, and J. M. Bowman
                                 "Quasiclassical Trajectory Calculations of the Rate Constant
                                 of the OH + HBr → Br + H2O Reaction Using a Full-Dimensional 
                                 AbInitio Potential Energy Surface Over the Temperature Range
                                 5 to 500 K",
                                 J. Phys. Chem. Lett. 5, 706-712 (2014).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_oxygen=0
    n_bromine=0
    natoms=len(geom)

    if natoms!=4:
        print("number of atoms not equal 4")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Br'):
                 n_bromine=n_bromine+1
            elif (geom[iatom][0]=='O'):
                 n_oxygen=n_oxygen+1
        if n_hydrogen!=2:
            print("number of Hydrogen atoms not equal 2 for H2OBr system")
            sys.exit()
        if n_bromine!=1:
            print("number of Bromine atoms not equal 1 for H2OBr system")
            sys.exit()
        if n_oxygen!=1:
            print("number of Oxygen atoms not equal 1 for H2OBr system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for H2OBr_PIP: input Cartesian should in order of HHOBr
    if surface=='H2OBr_PIP':
        geom_ordered=sorted(geom, key=lambda x: ("H", "O", "Br").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "O", "Br").index(x[0]))
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
