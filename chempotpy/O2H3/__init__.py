'''
O2H3
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for O2H3:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   O2H3_PIPNN:         single-state surface of O2H3, ground state,      P


        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for O2H3:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   O2H3_PIPNN:        single-state surface of O2H3, ground state,
                            availability: potential energy
                            functional form: PIP
                            corresponding surface in POTLIB: O2H3.zip
                            special emphsize on: OH + H2O → H2O + OH reaction 
                            ref: M. Bai, D. Lu, J. Li, 
                                 "Quasi-classical trajectory studies on the
                                 full-dimensional accurate potential energy 
                                 surface for the OH + H2O = H2O + OH reaction",
                                 Phys. Chem. Chem. Phys. 19, 17718-17725 (2017)


        """)

    return

def check(system, surface, geom):

    n_oxygen=0
    n_hydrogen=0
    natoms=len(geom)
    
    if natoms!=5:
        check=0
        print("number of atoms not equal 5")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='O'):
                 n_oxygen=n_oxygen+1
        if n_hydrogen!=3:
            print("number of Hydrogen atoms not equal 3 for O2H3 system")
            sys.exit()
        if n_oxygen!=2:
            print("number of Oxygen atoms not equal 2 for O2H3 system") 
            
    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for H3O2_BLL_2017: input Cartesian should in order of HHHOO
    if surface=='O2H3_PIPNN':
        geom_ordered=sorted(geom, key=lambda x: (x[0], x[1]))
        rank_ordered=sorted(rank, key=lambda x: (x[0], x[1]))
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
