'''
SiH5
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for SiH5:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   SiH5_LEP_1998v1:   single-state surface of SiH5, ground state,    P/G
    2.   SiH5_LEP_1998v2:   single-state surface of SiH5, ground state,    P/G

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for SiH5:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   SiH5_LEP_1998v1:   single-state surface of SiH5, ground state,
                            availability: potential energy
                            functional form: LEP
                            corresponding surface in POTLIB: sih5-1-2001.f
                            special emphsize on: SiH4 + H → SiH3 + H2 reaction 
                            ref: J. Espinosa-Garcia, J. Sanson, and J.C. Corchado,
                                 "The SiH4 + H → SiH3 + H2 reaction: Potential energy 
                                 surface, rate constants, and kinetic isotope 
                                 effects",
                                 J. Chem. Phys. 109, 466-473 (1998).
    ============================================================================================
    1.   SiH5_LEP_1998v2:   single-state surface of SiH5, ground state,
                            availability: potential energy
                            functional form: LEP
                            corresponding surface in POTLIB: sih5-2-2001.f
                            special emphsize on: SiH4 + H → SiH3 + H2 reaction 
                            ref: J. Espinosa-Garcia, J. Sanson, and J.C. Corchado,
                                 "The SiH4 + H → SiH3 + H2 reaction: Potential energy 
                                 surface, rate constants, and kinetic isotope 
                                 effects",
                                 J. Chem. Phys. 109, 466-473 (1998).



        """)

    return

def check(system, surface, geom):

    n_silicon=0
    n_hydrogen=0
    natoms=len(geom)

    if natoms!=6:
        print("number of atoms not equal 6")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Si'):
                 n_silicon=n_silicon+1
        if n_hydrogen!=5:
            print("number of Hydrogen atoms not equal 5 for SiH5 system")
            sys.exit()
        if n_silicon!=1:
            print("number of Silicon atoms not equal 1 for SiH5 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for SiH5_LEP_2000v1, SiH5_LEP_2000v2: input Cartesian should in order of HSiHHHH
    if surface=='SiH5_LEP_1998v1' or surface=='SiH5_LEP_1998v2':
        geom_ordered=sorted(geom, key=lambda x: ("Si", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Si", "H").index(x[0]))
        tmp_row=geom_ordered[0].copy()
        geom_ordered[0]=geom_ordered[1]
        geom_ordered[1]=tmp_row
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
