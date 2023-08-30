'''
ClNH3
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for ClNH3:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available. 

    Single-State Surfaces:
    1.  ClNH3_VBMM:           single-state surface of ClNH3, ground state,      P/G


        """)

    return

def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for ClNH3:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.  ClNH3_VBMM:           single-state surface of ClNH3, ground state,
                              availability: potential energy
                              functional form: VBMM
                              corresponding surface in POTLIB: clnh3-2012.f
                              special emphsize on: NH3 + Cl → NH2 + HCl reaction
                              ref: M. Monge-Palacios, C. Rangel,  J. C. Corchado,
                                   and J. Espinosa-Garcia, 
                                   "Analytical potential energy surface for the 
                                   reaction with intermediate complexes NH3 + Cl 
                                   → NH2 + HCl: Application to the kinetics study",
                                   Int. J. Quantum Chem. 112, 1887-1903 (2012).

        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_chlorine=0
    n_nitrogen=0
    natoms=len(geom)

    if natoms!=5:
        print("number of atoms not equal 5 for ClNH3 system")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='Cl'):
                 n_chlorine=n_chlorine+1
            elif (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='N'):
                 n_nitrogen=n_nitrogen+1
        if n_chlorine!=1:
            print("number of Chlorine atoms not equal 1 for ClNH3 system")
            sys.exit()
        if n_hydrogen!=3:
            print("number of Hydrogen atoms not equal 3 for ClNH3 system")
            sys.exit()
        if n_nitrogen!=1:
            print("number of Nitrogen atoms not equal 1 for ClNH3 system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for ClNH3_VBMM: input Cartesian should in order of HNHHCl
    if surface=='ClNH3_VBMM':
        geom_ordered=sorted(geom, key=lambda x: ("N", "H", "Cl").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("N", "H", "Cl").index(x[0]))
        for iatom in range(natoms):
            for idir in range(3):
                xyz[iatom][idir]=geom_ordered[iatom][idir+1]
        #now do re-ordering. 
        tmp_row=np.copy(xyz[0])
        xyz[0]=xyz[1]
        xyz[1]=tmp_row
        tmp_row=np.copy(rank_ordered[0])
        rank_ordered[0]=rank_ordered[1]
        rank_ordered[1]=tmp_row


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
