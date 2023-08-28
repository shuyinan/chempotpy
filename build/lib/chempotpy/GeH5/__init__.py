'''
GeH5
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for GeH5:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   GeH5_LEPS:         single-state surface of GeH5, ground state,      P/G

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule:
    The H atom furthest to Ge will be the first H, i.e. the H not in GeH4. 

        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for GeH5:
    LEPS: London–Eyring–Polanyi-Sato function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   GeH5_LEPS:         single-state surface of GeH5, ground state,
                            availability: potential energy
                            functional form: LEP
                            corresponding surface in POTLIB: geh5-2001.f
                            special emphsize on: GeH4 + H → GeH3 + H2 reaction
                            ref: J. Espinosa-Garcia,
                                 "Analytical potential energy surface for the 
                                 GeH4 + H → GeH3 + H2 reaction: Thermal and 
                                 vibrational-state selected rate constants and 
                                 kinetic isotope effects",
                                 J. Chem. Phys. 111, 9330-9336 (1999).

    ==VERY IMPORTANT NOTICE==
    the automatic ordering of input coordinates follows the folloiwng rule:
    The H atom furthest to Ge will be the first H, i.e. the H not in GeH4.

        """)

    return

def check(system, surface, geom):

    n_germanium=0
    n_hydrogen=0
    natoms=len(geom)

    if natoms!=6:
        print("number of atoms not equal 6")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Ge'):
                 n_germanium=n_germanium+1
        if n_hydrogen!=5:
            print("number of Hydrogen atoms not equal 5 for GeH5 system")
            sys.exit()
        if n_germanium!=1:
            print("number of Germanium atoms not equal 1 for GeH5 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for GeH5_LEPS: input Cartesian should in order of HGeHHHH
    if surface=='GeH5_LEPS':
        geom_ordered=sorted(geom, key=lambda x: ("Ge", "H").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("Ge", "H").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        r=np.zeros((5))
        for i in range(5):
            r[i]=np.linalg.norm(xyz[i+1]-xyz[0])
        max_idx=np.argmax(r)
        tmp_row=np.copy(xyz[1])
        xyz[1]=xyz[max_idx+1]
        xyz[max_idx+1]=tmp_row
        tmp_row=np.copy(rank_ordered[1])
        rank_ordered[1]=rank_ordered[max_idx+1]
        rank_ordered[max_idx+1]=tmp_row
        tmp_row=np.copy(xyz[1])
        xyz[1]=xyz[0]
        xyz[0]=tmp_row
        tmp_row=np.copy(rank_ordered[1])
        rank_ordered[1]=rank_ordered[0]  
        rank_ordered[0]=tmp_row 


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
