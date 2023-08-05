'''
H3Cl
'''
import numpy as np
import sys

def intro():

    print("""
    SHORT List of potential energy surfaces for H3Cl:
    P: potential energy, G: Gradient, D: Nonadiabatic coupling vector 
    This is listed at the end of description to show the availability.
    For example: P/G means both potential energy and gradient are available.

    Single-State Surfaces:
    1.   H3Cl_PIPNN:        single-state surface of H3Cl, ground state,      P


        """)

    return


def intro_detail():

    print("""
    DETAILED List of potential energy surfaces for H3Cl:
    LEP: London–Eyring–Polanyi function
    VBMM: valence bond molecular mechanics
    PIP: permutationally invariant polynomials
    PIP-NN: permutationally invariant polynomials followed by neural network 
    DDNN: diabatization by deep neural network 

    Single-State Surfaces:
    1.   H3Cl_PIPNN:        single-state surface of H3Cl, ground state,
                            availability: potential energy
                            functional form: PIP-NN
                            corresponding surface in POTLIB: H3Cl.zip
                            special emphsize on: H2 + HCl → H2 + HCl
                            ref: Q. Yao, M. Morita, C. Xie, N. Balakrishnan, and H. Guo
                                 "Globally Accurate Full-Dimensional Potential Energy 
                                 Surface for H2 + HCl Inelastic Scattering"
                                 J. Phys. Chem. A 123, 6578-6586 (2019).


        """)

    return

def check(system, surface, geom):

    n_hydrogen=0
    n_chlorine=0
    natoms=len(geom)

    if natoms!=4:
        print("number of atoms not equal 4")
        sys.exit()
    else:
        for iatom in range(natoms):
            if (geom[iatom][0]=='H'):
                 n_hydrogen=n_hydrogen+1
            elif (geom[iatom][0]=='Cl'):
                 n_chlorine=n_chlorine+1
        if n_hydrogen!=3:
            print("number of Hydrogen atoms not equal 3 for H3Cl system")
            sys.exit()
        if n_chlorine!=1:
            print("number of Chlorine atoms not equal 1 for H3Cl system")
            sys.exit()

    xyz=np.zeros((natoms,3))

    rank=np.copy(geom)
    for iatom in range(natoms):
        rank[iatom][1]=int(iatom)

    #for H3Cl_PIPNN: input Cartesian should in order of HHHCl
    if surface=='H3Cl_PIPNN':
        geom_ordered=sorted(geom, key=lambda x: ("H", "Cl").index(x[0]))
        rank_ordered=sorted(rank, key=lambda x: ("H", "Cl").index(x[0]))
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
